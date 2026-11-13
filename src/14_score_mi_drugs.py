# -*- coding: utf-8 -*-
"""14 第二层评估：用主模型给 MI 药打分，并在"MI 药 + 诱饵池"合集中算百分位与富集。

输入：data/processed/ligand_features_ecfp4.npz、_desc.csv、scaffold_splits.csv、
      results/tables/pool_decoys.csv、mi_gold_drugs.csv、ligand_model_perf.csv（主模型名）
输出：results/tables/mi_drug_scores.csv、logs/p7_experiments.md（追加）
运行：C:/Anaconda/Anaconda3/python.exe src/14_score_mi_drugs.py

MI 药可用数口径：161 -> 有结构 123 -> 标准化后可用 108（见 logs/pool_log.md），本脚本只用这 108 个。
结果无论正负照实报，禁止为改善它调任何参数。
"""
from pathlib import Path
import sys
import datetime

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdFingerprintGenerator
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import (DESC_COLS, LOGS, MAIN_TARGETS, SEED, TAB,  # noqa: E402
                          load_features, load_splits, make_model)
from pool_utils import prepare_smiles  # noqa: E402

RDLogger.DisableLog("rdApp.*")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def featurize(smiles_list):
    """返回 (成功行号列表, fps, descs, 失败记录)。身份一律用行号/ID 连接，不靠截断匹配。"""
    idx, fps, descs, failed = [], [], [], []
    for i, smi in enumerate(smiles_list):
        mol, reason = prepare_smiles(smi)
        if mol is None:
            failed.append({"row_index": i, "reason": reason})
            continue
        idx.append(i)
        fps.append(FPGEN.GetFingerprintAsNumPy(mol).astype(np.float32))
        descs.append([Descriptors.MolWt(mol), Crippen.MolLogP(mol),
                      rdMolDescriptors.CalcNumHBD(mol), rdMolDescriptors.CalcNumHBA(mol),
                      rdMolDescriptors.CalcTPSA(mol), rdMolDescriptors.CalcNumRotatableBonds(mol),
                      rdMolDescriptors.CalcNumRings(mol), rdMolDescriptors.CalcNumAromaticRings(mol),
                      rdMolDescriptors.CalcFractionCSP3(mol)])
    if not idx:
        return np.array([], dtype=int), np.zeros((0, 2048), np.float32), np.zeros((0, 9)), failed
    return np.array(idx), np.vstack(fps), np.array(descs, dtype=np.float32), failed


def main():
    ecfp, desc = load_features()
    splits = load_splits()
    perf = pd.read_csv(TAB / "ligand_model_perf.csv", encoding="utf-8-sig")
    main_model = perf.loc[perf["is_main"], "model"].iloc[0]

    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")
    gold = gold[gold["smiles"].notna() & (gold["smiles"].astype(str).str.strip() != "")]
    g_idx, G_fp, G_desc, g_failed = featurize(gold["smiles"].tolist())
    # 用行号显式连接，保证分子身份与预测分数一一对应
    gold = gold.iloc[g_idx].reset_index(drop=True)
    pd.DataFrame(g_failed).to_csv(TAB / "mi_drug_featurize_failed.csv", index=False,
                                  encoding="utf-8-sig")
    print(f"MI 药物可用（标准化后）：{len(gold)}；特征化失败 {len(g_failed)}（已记录）")

    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")

    out_rows, log_rows = [], []
    for target in MAIN_TARGETS:
        df = splits[splits["target"] == target].reset_index(drop=True)
        ids = df["molecule_chembl_id"].tolist()
        y = df["y"].to_numpy()
        X = np.vstack([ecfp[m] for m in ids]).astype(np.float32)
        model = make_model(main_model, y)
        model.fit(X, y)

        d_ids = dec.loc[dec["target"] == target, "molecule_chembl_id"].tolist()
        Xd = np.vstack([ecfp[m] for m in d_ids]).astype(np.float32)
        p_dec = model.predict_proba(Xd)[:, 1]
        p_gold = model.predict_proba(G_fp)[:, 1]

        allp = np.concatenate([p_gold, p_dec])
        pct = np.array([(allp <= v).mean() * 100 for v in p_gold])

        u, pval = mannwhitneyu(p_gold, p_dec, alternative="two-sided")
        k = max(1, int(round(0.01 * len(allp))))
        top_idx = np.argsort(-allp)[:k]
        n_gold_top = int((top_idx < len(p_gold)).sum())
        enrich = (n_gold_top / k) / (len(p_gold) / len(allp))

        for i, r in gold.iterrows():
            out_rows.append({"target": target, "molecule_chembl_id": r["molecule_chembl_id"],
                             "pref_name": r["pref_name"], "prob": float(p_gold[i]),
                             "percentile_rank": float(pct[i])})

        log_rows.append((target, float(np.median(p_gold)), float(np.median(p_dec)),
                         float(u), float(pval), enrich, n_gold_top, k))
        print(f"{target}: MI 中位 {np.median(p_gold):.3f} vs 诱饵中位 {np.median(p_dec):.3f}，"
              f"p={pval:.3e}，top1% 富集 {enrich:.2f}x")

    out = pd.DataFrame(out_rows)
    out.to_csv(TAB / "mi_drug_scores.csv", index=False, encoding="utf-8-sig")

    with open(LOGS / "p7_experiments.md", "a", encoding="utf-8") as f:
        f.write(f"\n## {datetime.datetime.now().isoformat(timespec='seconds')} MI 药评估（14）\n")
        f.write(f"- 主模型：{main_model}；MI 药可用 {len(gold)} 个（161 -> 有结构 123 -> 可用 {len(gold)}）\n")
        f.write("- 方向说明：U 统计量偏大表示 MI 药预测概率高于诱饵；enrich>1 表示 top1% 富集\n")
        for t, mg, md, u, p, e, ng, k in log_rows:
            f.write(f"- {t}: MI 中位 {mg:.4f} vs 诱饵中位 {md:.4f}，"
                    f"Mann-Whitney U={u:.1f}, p={p:.3e}，top1% 富集 {e:.2f}x（{ng}/{k}）\n")
        f.write("")

    print("[14_score_mi_drugs] 完成")


if __name__ == "__main__":
    main()