# -*- coding: utf-8 -*-
"""29 最终图与统计：库得分统计、图8（含控制药物标注）、图9（逐靶点排除后重绘）。

输出：results/tables/library_score_stats.csv、reports/fig8_score_property.png、
      reports/fig9_mi_enrichment.png
运行：C:/Anaconda/Anaconda3/python.exe src/29_final_figures.py
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import MAIN_TARGETS, ROOT, TAB  # noqa: E402

REP = ROOT / "reports"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
                     "figure.dpi": 100, "savefig.dpi": 300})

CONTROLS = {"CRISABOROLE": "CHEMBL484785", "IBUDILAST": "CHEMBL19449",
            "TAVABOROLE": "CHEMBL443052"}


def train_keys(target):
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    return (set(act.loc[act["target"] == target, "inchikey14"])
            | set(dec.loc[dec["target"] == target, "inchikey14"]))


def main():
    sc = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    pool = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")

    seen = set()
    for t in MAIN_TARGETS:
        seen |= train_keys(t)
    cand = sc[(sc["is_mi_gold"] == 0) & (~sc["inchikey14"].isin(seen))].copy()
    cand = cand.merge(pool[["molecule_chembl_id", "mw"]], on="molecule_chembl_id", how="left")

    # ---- 5.7 分布统计（numpy 线性插值分位数，method='linear'）
    p = cand["prob_PDE4B"].to_numpy()
    stats = pd.DataFrame([{
        "n_candidates": len(p),
        "median": float(np.median(p)),
        "p99_method_linear": float(np.percentile(p, 99, method="linear")),
        "max": float(p.max()),
        "n_above_0_5": int((p > 0.5).sum()),
    }])
    stats.to_csv(TAB / "library_score_stats.csv", index=False, encoding="utf-8-sig")
    print(stats.to_string(index=False))

    # ---- 图 8：分布 + 得分-分子量，标注控制药物
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    axes[0].hist(p, bins=60, color="#4c72b0")
    axes[0].axvline(stats["median"].iloc[0], ls="--", lw=1, color="#333333",
                    label=f"median = {stats['median'].iloc[0]:.3f}")
    axes[0].axvline(stats["p99_method_linear"].iloc[0], ls=":", lw=1, color="#c44e52",
                    label=f"99th percentile = {stats['p99_method_linear'].iloc[0]:.3f}")
    axes[0].set_xlabel("PDE4B model score")
    axes[0].set_ylabel("number of molecules")
    axes[0].set_title(f"(A) Score distribution (n = {len(p)})")
    axes[0].legend(frameon=False)

    axes[1].scatter(cand["mw"], cand["prob_PDE4B"], s=7, alpha=0.35, color="#4c72b0")
    for name, cid in CONTROLS.items():
        row = cand[cand["molecule_chembl_id"] == cid]
        if row.empty:
            continue
        x, y = float(row["mw"].iloc[0]), float(row["prob_PDE4B"].iloc[0])
        axes[1].scatter([x], [y], s=55, facecolors="none", edgecolors="#d62728", linewidths=1.4)
        axes[1].annotate(name, (x, y), textcoords="offset points", xytext=(6, 6), fontsize=8,
                         color="#d62728")
    axes[1].set_xlabel("molecular weight (Da)")
    axes[1].set_ylabel("PDE4B model score")
    axes[1].set_title("(B) Model score vs molecular weight")
    fig.tight_layout()
    fig.savefig(REP / "fig8_score_property.png", dpi=300)
    plt.close(fig)

    # ---- 图 9：逐靶点排除后重绘（与 22 号规则一致）
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")
    from pool_utils import descriptors, prepare_smiles
    gold_ik = set()
    for s in gold["smiles"].dropna():
        mol, _ = prepare_smiles(s)
        if mol is not None:
            gold_ik.add(descriptors(mol)["inchikey14"])
    sc["is_mi"] = sc["inchikey14"].isin(gold_ik).astype(int)

    panels = [(t, f"prob_{t}", [t]) for t in MAIN_TARGETS] + \
             [("consensus", "prob_mean", MAIN_TARGETS)]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, (name, col, tl) in zip(axes.ravel(), panels):
        keys = set()
        for t in tl:
            keys |= train_keys(t)
        sub = sc[~sc["inchikey14"].isin(keys)]
        mi = sub.loc[sub.is_mi == 1, col].to_numpy()
        non = sub.loc[sub.is_mi == 0, col].to_numpy()
        bp = ax.boxplot([non, mi], tick_labels=[f"non-MI (n={len(non)})", f"MI (n={len(mi)})"],
                        showfliers=False, patch_artist=True, widths=0.5)
        for patch, c in zip(bp["boxes"], ["#9ecae1", "#f4a582"]):
            patch.set_facecolor(c)
        for med in bp["medians"]:
            med.set_color("#333333")
        ax.set_title(name)
        ax.set_ylabel("model score")
        ax.grid(axis="y", ls=":", alpha=0.4)
    fig.suptitle("Model scores: MI-indication molecules vs other approved drugs "
                 "(training overlap excluded per target)")
    fig.tight_layout()
    fig.savefig(REP / "fig9_mi_enrichment.png", dpi=300)
    plt.close(fig)
    print("[29] 图 8/图 9 与库统计已生成")


if __name__ == "__main__":
    main()