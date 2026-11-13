# -*- coding: utf-8 -*-
"""28 口径核对与补充分析重生成：候选集合、共识排名（mean percentile 为主，Borda 为辅）、
控制药物排名（PDE4B 单靶点在共同候选集合内）、图 8 双面板。

输出：reports/supp_S1_rank_consensus_top25.csv、supp_S1_rank_consensus_methods.csv、
      supp_S2_score_property_correlation.csv、supp_S3_control_molecules.csv、
      reports/fig8_score_property.png、logs/scope_check.md
运行：C:/Anaconda/Anaconda3/python.exe src/28_supplements_v2.py
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import MAIN_TARGETS, ROOT, TAB, LOGS  # noqa: E402

REP = ROOT / "reports"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "figure.dpi": 100, "savefig.dpi": 300})


def main():
    sc = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    pool = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    # 统一规则：共识候选只排除三个参与评分模型（PDE4B/PFKFB3/EPAS1）的训练数据；
    # PDK4 仅作探索性基准，不参与共识排除
    seen = set(act.loc[act["target"].isin(MAIN_TARGETS), "inchikey14"]) | \
        set(dec.loc[dec["target"].isin(MAIN_TARGETS), "inchikey14"])

    n0 = len(sc)
    n_mi = int((sc["is_mi_gold"] == 1).sum())
    non_mi = sc[sc["is_mi_gold"] == 0]
    n_overlap = int(non_mi["inchikey14"].isin(seen).sum())
    cand = non_mi[~non_mi["inchikey14"].isin(seen)].copy()
    print(f"打分表 {n0} | MI 标记 {n_mi} | 非 MI {len(non_mi)} | 训练重叠 {n_overlap} "
          f"| 候选 {len(cand)}")

    # 与短名单脚本口径逐行比对
    sl = pd.read_csv(TAB / "repurposing_shortlist.csv", encoding="utf-8-sig")
    sl_ids = set(sl["chembl_id"])
    cand_ids = set(cand["molecule_chembl_id"])
    only_sl = sl_ids - cand_ids
    print(f"短名单候选差异：仅短名单 {len(only_sl)} 个 {sorted(only_sl)[:5]}")

    # 性质缺口
    m = cand.merge(pool[["molecule_chembl_id", "max_phase", "mw", "clogp"]],
                   on="molecule_chembl_id", how="left")
    n_no_prop = int(m["mw"].isna().sum())
    print(f"候选合并性质后缺失 mw 的分子：{n_no_prop}")
    if n_no_prop:
        print(m.loc[m["mw"].isna(), ["molecule_chembl_id", "pref_name"]].to_string(index=False))

    # ---- 共识排名：mean percentile（主要）与 Borda（补充）
    m["pct_PDE4B"] = m["prob_PDE4B"].rank(pct=True)
    m["pct_PFKFB3"] = m["prob_PFKFB3"].rank(pct=True)
    m["pct_EPAS1"] = m["prob_EPAS1"].rank(pct=True)
    m["mean_percentile"] = m[["pct_PDE4B", "pct_PFKFB3", "pct_EPAS1"]].mean(axis=1)
    m["rank_mean_pct"] = m["mean_percentile"].rank(ascending=False, method="min")
    m["rank_prob_mean"] = m["prob_mean"].rank(ascending=False, method="min")
    for t in MAIN_TARGETS:
        m[f"rank_{t}"] = m[f"prob_{t}"].rank(ascending=False, method="min")
    m["rank_sum"] = m[[f"rank_{t}" for t in MAIN_TARGETS]].sum(axis=1)
    # Borda：秩和越小越好，故按升序排名（此前实现按降序，方向相反，已修正）
    m["rank_borda"] = m["rank_sum"].rank(ascending=True, method="min")

    base20 = set(m.nsmallest(20, "rank_mean_pct")["molecule_chembl_id"])
    rho_pm = spearmanr(m["rank_mean_pct"], m["rank_prob_mean"]).correlation
    rho_bd = spearmanr(m["rank_mean_pct"], m["rank_borda"]).correlation
    methods = pd.DataFrame([
        {"comparison": "mean_percentile vs prob_mean", "spearman_rho": round(float(rho_pm), 4),
         "top20_overlap": len(base20 & set(m.nsmallest(20, "rank_prob_mean")["molecule_chembl_id"])),
         "top50_overlap": len(set(m.nsmallest(50, "rank_mean_pct")["molecule_chembl_id"])
                              & set(m.nsmallest(50, "rank_prob_mean")["molecule_chembl_id"]))},
        {"comparison": "mean_percentile vs borda_rank_sum", "spearman_rho": round(float(rho_bd), 4),
         "top20_overlap": len(base20 & set(m.nsmallest(20, "rank_borda")["molecule_chembl_id"])),
         "top50_overlap": len(set(m.nsmallest(50, "rank_mean_pct")["molecule_chembl_id"])
                              & set(m.nsmallest(50, "rank_borda")["molecule_chembl_id"]))},
    ])
    methods.to_csv(REP / "supp_S1_rank_consensus_methods.csv", index=False, encoding="utf-8-sig")
    print(methods.to_string(index=False))

    top25 = m.sort_values("rank_mean_pct").head(25)
    out = top25[["pref_name", "molecule_chembl_id", "max_phase", "prob_PDE4B", "prob_PFKFB3",
                 "prob_EPAS1", "prob_mean", "mean_percentile", "rank_mean_pct",
                 "rank_prob_mean", "rank_borda"]].copy()
    out.insert(0, "rank", range(1, len(out) + 1))
    out.to_csv(REP / "supp_S1_rank_consensus_top25.csv", index=False, encoding="utf-8-sig")

    # ---- 性质相关
    rows = []
    for t in MAIN_TARGETS + ["consensus"]:
        col = "prob_mean" if t == "consensus" else f"prob_{t}"
        for prop in ["mw", "clogp"]:
            d = m[[col, prop]].dropna()
            rho, p = spearmanr(d[col], d[prop])
            rows.append({"score": t, "property": prop, "spearman_rho": round(float(rho), 4),
                         "p_value": float(p), "n": len(d)})
    s2 = pd.DataFrame(rows)
    s2.to_csv(REP / "supp_S2_score_property_correlation.csv", index=False, encoding="utf-8-sig")

    # ---- 控制药物：统一在共同候选集合内算 PDE4B 排名
    controls = {
        "CRISABOROLE": ("CHEMBL484785", "PDE4 inhibitor (approved, topical)"),
        "ROFLUMILAST": ("CHEMBL193240", "PDE4 inhibitor (approved)"),
        "APREMILAST": ("CHEMBL514800", "PDE4 inhibitor (approved)"),
        "IBUDILAST": ("CHEMBL19449", "PDE4 inhibitor (approved, Japan)"),
        "TAVABOROLE": ("CHEMBL443052", "LeuRS inhibitor (antifungal), not a PDE4 inhibitor"),
        "CILOSTAZOL": ("CHEMBL799", "PDE3 inhibitor (negative control for PDE4)"),
    }
    train_pos_ids = set(act["molecule_chembl_id"])
    rows3 = []
    for name, (cid, note) in controls.items():
        in_cand = m[m["molecule_chembl_id"] == cid]
        if in_cand.empty:
            in_train = cid in train_pos_ids
            rows3.append({"pref_name": name, "chembl_id": cid, "annotation": note,
                          "inchikey14": "", "status": ("excluded: training positive" if in_train
                                                        else "absent from scoring table"),
                          "in_training": bool(in_train), "prob_PDE4B": np.nan,
                          "rank_PDE4B": np.nan, "rank_PDE4B_denom": np.nan,
                          "tied_with": np.nan, "rank_consensus_mean_pct": np.nan})
            continue
        r = in_cand.iloc[0]
        ties = int((m["prob_PDE4B"] == r["prob_PDE4B"]).sum())
        rows3.append({"pref_name": name, "chembl_id": cid, "annotation": note,
                      "inchikey14": r["inchikey14"], "status": "independent (not in training)",
                      "in_training": bool(r["inchikey14"] in seen),
                      "prob_PDE4B": round(float(r["prob_PDE4B"]), 4),
                      "rank_PDE4B": int(r["rank_PDE4B"]),
                      "rank_PDE4B_denom": len(m), "tied_with": ties,
                      "rank_consensus_mean_pct": int(r["rank_mean_pct"])})
    s3 = pd.DataFrame(rows3)
    s3.to_csv(REP / "supp_S3_control_molecules.csv", index=False, encoding="utf-8-sig")
    print(s3[["pref_name", "status", "prob_PDE4B", "rank_PDE4B", "rank_PDE4B_denom",
              "rank_consensus_mean_pct"]].to_string(index=False))

    # ---- 图 8：同一候选集合的双面板
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    axes[0].hist(m["prob_PDE4B"], bins=50, color="#4c72b0")
    axes[0].set_xlabel("PDE4B model score")
    axes[0].set_ylabel("number of molecules")
    axes[0].set_title(f"Score distribution (n = {len(m)})")
    axes[1].scatter(m["mw"], m["prob_PDE4B"], s=8, alpha=0.45, color="#4c72b0")
    axes[1].set_xlabel("molecular weight (Da)")
    axes[1].set_ylabel("PDE4B model score")
    axes[1].set_title(f"Score vs molecular weight (Spearman rho = "
                      f"{s2[(s2.score=='PDE4B')&(s2.property=='mw')]['spearman_rho'].iloc[0]:.3f})")
    fig.tight_layout()
    fig.savefig(REP / "fig8_score_property.png", dpi=300)
    plt.close(fig)

    (LOGS / "scope_check.md").write_text(
        f"""# 口径核对

- 打分表 {n0}；MI 标记 {n_mi}；非 MI {len(non_mi)}；与训练数据重叠 {n_overlap}；**共同候选 {len(m)}**。
- 与 `repurposing_shortlist.csv` 的候选差集：{len(only_sl)} 个。
- 合并性质后缺失 mw 的分子：{n_no_prop}（缺者即为图注需说明的分子）。
- 排名定义：PDE4B 排名在共同候选集合内按 PDE4B 得分降序，并列取同分同名次（method="min"），分母 {len(m)}。
- 聚合：主要用 mean percentile（三靶点百分位均值），Borda（秩和）仅作补充比较。

{methods.to_string(index=False)}
""", encoding="utf-8")
    print("[28] 完成")


if __name__ == "__main__":
    main()