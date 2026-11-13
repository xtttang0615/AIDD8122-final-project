# -*- coding: utf-8 -*-
"""24 P8 图：MI 药物富集分布 + 短名单热图（300 dpi，字号>=8pt）。

输入：results/tables/approved_scores.csv、repurposing_shortlist.csv
输出：results/figures/fig7_mi_enrichment.png、fig7b_shortlist_heatmap.png
运行：C:/Anaconda/Anaconda3/python.exe src/24_p8_figures.py
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

FIG = ROOT / "results" / "figures"
plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 100, "savefig.dpi": 300,
})


def plot_enrichment(sc):
    panels = [(t, f"prob_{t}") for t in MAIN_TARGETS] + [("consensus", "prob_mean")]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    for ax, (name, col) in zip(axes.ravel(), panels):
        mi = sc.loc[sc.is_mi_gold == 1, col].dropna().to_numpy()
        non = sc.loc[sc.is_mi_gold == 0, col].dropna().to_numpy()
        bp = ax.boxplot([non, mi], labels=[f"non-MI (n={len(non)})", f"MI (n={len(mi)})"],
                        showfliers=False, patch_artist=True, widths=0.5)
        for patch, c in zip(bp["boxes"], ["#9ecae1", "#f4a582"]):
            patch.set_facecolor(c)
        for med in bp["medians"]:
            med.set_color("#333333")
        ax.set_title(name)
        ax.set_ylabel("predicted probability")
        ax.grid(axis="y", ls=":", alpha=0.4)
    fig.suptitle("MI drugs vs other approved drugs (model scores)")
    fig.tight_layout()
    fig.savefig(FIG / "fig7_mi_enrichment.png", dpi=300)
    plt.close(fig)
    print("[fig] fig7_mi_enrichment.png")


def plot_heatmap(sl):
    top = sl.head(20)
    M = top[["prob_PDE4B", "prob_PFKFB3", "prob_EPAS1"]].to_numpy()
    labels = [str(n)[:22] for n in top["pref_name"]]
    fig, ax = plt.subplots(figsize=(6.5, 8))
    im = ax.imshow(M, aspect="auto", cmap="YlOrRd")
    ax.set_xticks(range(3))
    ax.set_xticklabels(["PDE4B", "PFKFB3", "EPAS1"], rotation=20)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=8)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            ax.text(j, i, f"{M[i, j]:.3f}", ha="center", va="center", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("predicted probability")
    ax.set_title("Top 20 repurposing shortlist")
    fig.tight_layout()
    fig.savefig(FIG / "fig7b_shortlist_heatmap.png", dpi=300)
    plt.close(fig)
    print("[fig] fig7b_shortlist_heatmap.png")


def main():
    sc = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    sl = pd.read_csv(TAB / "repurposing_shortlist.csv", encoding="utf-8-sig")
    # MI 富集图（fig9）由 src/29_final_figures.py 统一生成（按靶点排除训练重叠），此处不再重复生成
    plot_heatmap(sl)
    print("[24_p8_figures] 完成（仅短名单热图；fig9 见 29_final_figures.py）")


if __name__ == "__main__":
    main()