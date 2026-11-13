# -*- coding: utf-8 -*-
"""22 MI 药物富集检验（定稿版）：108 个 MI 金标准 vs 已上市药物池中的其余药物。

输入：results/tables/approved_scores.csv、mi_gold_drugs.csv
输出：results/tables/mi_enrichment_final.csv、logs/p8_log.md（追加）
运行：C:/Anaconda/Anaconda3/python.exe src/22_mi_enrichment.py

指标：Mann-Whitney U（单尾，MI 更高）、AUC-like = P(score_MI > score_nonMI)（并列计 0.5）、
top1%/top5% 富集倍数。结论照实写，p>0.05 时明确写"未观察到富集"。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, MAIN_TARGETS, TAB  # noqa: E402
from pool_utils import descriptors, prepare_smiles  # noqa: E402


def ef_at(y, p, frac):
    n = len(y)
    k = max(1, int(round(frac * n)))
    idx = np.argsort(-p)[:k]
    base = y.mean()
    return float((y[idx].sum() / k) / base) if base > 0 else np.nan


def train_keys(target):
    """该靶点模型的全部训练结构键（正例 pool_actives + 负例 pool_decoys）。"""
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    return (set(act.loc[act.target == target, "inchikey14"])
            | set(dec.loc[dec.target == target, "inchikey14"]))


def main():
    scores = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")

    gold_ik = set()
    for s in gold["smiles"].dropna():
        mol, _ = prepare_smiles(s)
        if mol is not None:
            gold_ik.add(descriptors(mol)["inchikey14"])
    print(f"MI 适应证分子（ChEMBL 记录，标准化去重后）：{len(gold_ik)}")

    scores["is_mi"] = scores["inchikey14"].isin(gold_ik).astype(int)

    rows = []
    cols = [(t, f"prob_{t}", {t}) for t in MAIN_TARGETS]
    cols.append(("consensus", "prob_mean", set(MAIN_TARGETS)))

    excluded_log = []
    for name, col, tl in cols:
        keys = set()
        for t in tl:
            keys |= train_keys(t)
        mask = ~scores["inchikey14"].isin(keys)
        sub = scores.loc[mask]
        n_excl = int((~mask).sum())
        y = sub["is_mi"].to_numpy()
        p = sub[col].to_numpy()
        mi, non = p[y == 1], p[y == 0]
        excluded_log.append(f"- {name}: 剔除与该模型训练数据重叠 {n_excl} 个分子，"
                            f"剩余评价集 {len(sub)}（MI {len(mi)} / 非 MI {len(non)}）")
        if len(mi) == 0 or len(non) == 0:
            rows.append({"target": name, "n_MI": len(mi), "n_nonMI": len(non),
                         "median_MI": np.nan, "median_nonMI": np.nan, "U_p": np.nan,
                         "auc_like": np.nan, "ef1": np.nan, "ef5": np.nan,
                         "conclusion": "评价集不足"})
            continue
        u, pval = mannwhitneyu(mi, non, alternative="greater")
        rows.append({"target": name, "n_MI": len(mi), "n_nonMI": len(non),
                     "median_MI": float(np.median(mi)), "median_nonMI": float(np.median(non)),
                     "U_p": float(pval), "auc_like": float(u / (len(mi) * len(non))),
                     "ef1": ef_at(y, p, 0.01), "ef5": ef_at(y, p, 0.05),
                     "conclusion": "观察到富集" if pval < 0.05 else "未观察到富集"})

    res = pd.DataFrame(rows)
    res.to_csv(TAB / "mi_enrichment_final.csv", index=False, encoding="utf-8-sig")
    print(res.round(4).to_string(index=False))

    with open(LOGS / "p8_log.md", "a", encoding="utf-8") as f:
        f.write("## 三（修订）、MI 适应证关联检查（排除训练数据重叠后）\n\n")
        f.write("- 集合定义：ChEMBL 中具有 MI 适应证记录的分子（非“已获批 AMI 药物金标准”）\n")
        f.write("\n".join(excluded_log) + "\n\n")
        f.write(res.round(4).to_string(index=False) + "\n\n")
    print("[22_mi_enrichment] 完成")


if __name__ == "__main__":
    main()