# -*- coding: utf-8 -*-
"""21 用三个主靶点模型给已上市药物池打分（RF + ECFP4，参数与 P7 完全一致）。

输入：results/tables/approved_pool.csv、data/processed/scaffold_splits.csv、ligand_features_ecfp4.npz
输出：results/tables/approved_scores.csv、logs/p8_log.md（追加覆盖度）
运行：C:/Anaconda/Anaconda3/python.exe src/21_score_approved.py

不调参、不换特征；按 P7 日志中的参数重建 RF（random_state=42，结果与已保存模型一致）。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from rdkit import RDLogger
from rdkit.Chem import rdFingerprintGenerator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, MAIN_TARGETS, SEED, TAB, load_features, load_splits, make_model  # noqa: E402
from pool_utils import descriptors, prepare_smiles  # noqa: E402

RDLogger.DisableLog("rdApp.*")
FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def main():
    ecfp, _ = load_features()
    splits = load_splits()
    pool = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")

    # MI 金标准被 20 号脚本从候选池剔除，此处并入打分表（打标 is_mi_gold），供 22 号做富集比较
    gold_recs = []
    for r in gold.itertuples(index=False):
        if not isinstance(r.smiles, str) or not r.smiles.strip():
            continue
        mol, _ = prepare_smiles(r.smiles)
        if mol is None:
            continue
        gold_recs.append({"molecule_chembl_id": r.molecule_chembl_id,
                          "pref_name": r.pref_name, "smiles": r.smiles,
                          "inchikey14": descriptors(mol)["inchikey14"]})
    gold_df = pd.DataFrame(gold_recs)
    gold_ik = set(gold_df["inchikey14"])

    allm = pd.concat([pool[["molecule_chembl_id", "pref_name", "smiles", "inchikey14"]],
                      gold_df], ignore_index=True)
    allm = allm.drop_duplicates("inchikey14", keep="first").reset_index(drop=True)
    allm["is_mi_gold"] = allm["inchikey14"].isin(gold_ik).astype(int)
    print(f"候选池 {len(pool)} + MI 金标准 {(allm['is_mi_gold'] == 1).sum()}（去重后共 {len(allm)}）")

    fps = []
    for s in allm["smiles"]:
        mol, _ = prepare_smiles(s)
        fps.append(FPGEN.GetFingerprintAsNumPy(mol).astype(np.float32))
    Xp = np.vstack(fps)
    print(f"打分矩阵：{Xp.shape}")

    out = allm[["molecule_chembl_id", "pref_name", "inchikey14", "is_mi_gold"]].copy()
    for target in MAIN_TARGETS:
        df = splits[splits["target"] == target].reset_index(drop=True)
        y = df["y"].to_numpy()
        X = np.vstack([ecfp[m] for m in df["molecule_chembl_id"]]).astype(np.float32)
        model = make_model("rf", y)
        model.fit(X, y)
        p = model.predict_proba(Xp)[:, 1]
        out[f"prob_{target}"] = p
        out[f"rank_{target}"] = pd.Series(-p).rank(method="min").astype(int)

    probs = out[[f"prob_{t}" for t in MAIN_TARGETS]]
    out["prob_mean"] = probs.mean(axis=1)
    out["prob_min_rank"] = out[[f"rank_{t}" for t in MAIN_TARGETS]].min(axis=1)
    out = out[["molecule_chembl_id", "pref_name", "inchikey14", "is_mi_gold",
               "prob_PDE4B", "prob_PFKFB3", "prob_EPAS1", "prob_mean", "prob_min_rank"]]
    out.to_csv(TAB / "approved_scores.csv", index=False, encoding="utf-8-sig")

    n_all = int(out[[f"prob_{t}" for t in MAIN_TARGETS]].notna().all(axis=1).sum())
    print(f"三靶点均可打分：{n_all}/{len(out)}（{n_all / len(out):.1%}）")

    with open(LOGS / "p8_log.md", "a", encoding="utf-8") as f:
        f.write(f"- 打分覆盖：{n_all}/{len(out)} 个分子三靶点均可打分（{n_all / len(out):.1%}）\n")
        f.write(f"- 打分表含候选池 {len(pool)} 个 + MI 金标准 "
                f"{int((out['is_mi_gold'] == 1).sum())} 个（后者仅供富集检验与短名单排除）\n")
        f.write("- 打分模型：RF(500)+ECFP4，与 P7 主模型参数一致，未调参\n\n")
    print("[21_score_approved] 完成")


if __name__ == "__main__":
    main()