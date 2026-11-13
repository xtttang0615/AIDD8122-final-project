# -*- coding: utf-8 -*-
"""30 报告表格核查：表7 双列表 top10、聚合列表重叠（含比较对象与交集 ID）。

输出：reports/table7_check.csv、reports/top_overlap_check.csv
运行：C:/Anaconda/Anaconda3/python.exe src/30_verify_tables.py
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import MAIN_TARGETS, ROOT, TAB  # noqa: E402

REP = ROOT / "reports"
TOPN = 20


def train_keys(t):
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    return (set(act.loc[act["target"] == t, "inchikey14"])
            | set(dec.loc[dec["target"] == t, "inchikey14"]))


def main():
    sc = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    pool = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")
    seen = set()
    for t in MAIN_TARGETS:
        seen |= train_keys(t)
    c = sc[(sc["is_mi_gold"] == 0) & (~sc["inchikey14"].isin(seen))].copy()
    c = c.merge(pool[["molecule_chembl_id", "max_phase"]], on="molecule_chembl_id", how="left")

    # 排名定义：降序，同分同名次（method="min"）
    c["rank_prob_mean"] = c["prob_mean"].rank(ascending=False, method="min")
    c["rank_PDE4B"] = c["prob_PDE4B"].rank(ascending=False, method="min")
    pct = c[[f"prob_{t}" for t in MAIN_TARGETS]].rank(pct=True)
    c["mean_percentile"] = pct.mean(axis=1)
    c["rank_mean_pct"] = c["mean_percentile"].rank(ascending=False, method="min")
    rs = c[[f"prob_{t}" for t in MAIN_TARGETS]].rank(ascending=False, method="min").sum(axis=1)
    c["borda"] = rs.rank(ascending=True, method="min")

    # 表7：两个 top10（固定取 10 个，按排名再按得分去重）
    a = c.nsmallest(10, "rank_prob_mean")[["rank_prob_mean", "pref_name", "molecule_chembl_id",
                                           "prob_PDE4B", "prob_PFKFB3", "prob_EPAS1", "prob_mean"]]
    b = c.nsmallest(10, "rank_PDE4B")[["rank_PDE4B", "pref_name", "molecule_chembl_id",
                                       "prob_PDE4B", "prob_mean"]]
    a.insert(0, "list", "mean_score")
    b.insert(0, "list", "PDE4B_only")
    t7 = pd.concat([a.rename(columns={"rank_prob_mean": "rank"}),
                    b.rename(columns={"rank_PDE4B": "rank"})], ignore_index=True)
    t7.to_csv(REP / "table7_check.csv", index=False, encoding="utf-8-sig")
    print("=== 表7 双列表 top10（n=1,823）===")
    print(t7.to_string(index=False))

    # 重叠：固定取前 N 个（同分不扩展），并报告边界同分
    def ids(col):
        return list(c.nsmallest(TOPN, col)["molecule_chembl_id"])

    rows = []
    for col_a, col_b, name in [
        ("rank_prob_mean", "rank_mean_pct", "mean_score vs mean_percentile"),
        ("rank_mean_pct", "rank_PDE4B", "mean_percentile vs PDE4B_only"),
        ("rank_prob_mean", "rank_PDE4B", "mean_score vs PDE4B_only"),
        ("rank_mean_pct", "borda", "mean_percentile vs Borda"),
    ]:
        A, B = set(ids(col_a)), set(ids(col_b))
        rows.append({"comparison": name, "set": "fixed_top20", "n_candidates": len(c),
                     "overlap": len(A & B), "overlap_ids": ";".join(sorted(A & B))})
        rows.append({"comparison": name, "set": "fixed_top50", "n_candidates": len(c),
                     "overlap": len(set(c.nsmallest(50, col_a)["molecule_chembl_id"])
                                    & set(c.nsmallest(50, col_b)["molecule_chembl_id"])),
                     "overlap_ids": ""})
    ov = pd.DataFrame(rows)
    ov.to_csv(REP / "top_overlap_check.csv", index=False, encoding="utf-8-sig")
    print("\n=== 重叠（固定取 N，同分不扩展）===")
    print(ov.to_string(index=False))


if __name__ == "__main__":
    main()