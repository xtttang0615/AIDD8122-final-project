# -*- coding: utf-8 -*-
"""质控图：PCA、每样本箱线图、密度图（300 dpi，字号>=8pt）。

输入：data/processed/{GSE}_gene.csv、data/processed/{GSE}_meta.csv
输出：results/figures/fig2a_pca_gse66360_by_role.png
      results/figures/fig2b_pca_gse66360_by_group.png
      results/figures/fig2c_pca_gse48060_by_group.png
      results/figures/fig2d_boxplot_gse66360.png
      results/figures/fig2e_density_gse66360.png
运行：C:/Anaconda/Anaconda3/python.exe src/04_qc_plot.py
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
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 100,
    "savefig.dpi": 300,
})


def load(gse):
    expr = pd.read_csv(PROCESSED / f"{gse}_gene.csv", index_col=0, encoding="utf-8-sig")
    meta = pd.read_csv(PROCESSED / f"{gse}_meta.csv", encoding="utf-8-sig").set_index("sample_id")
    meta = meta.loc[expr.columns]
    return expr, meta


def pca_plot(expr, meta, color_col, title, out, palette):
    X = expr.T.values
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)
    var = pca.explained_variance_ratio_ * 100
    fig, ax = plt.subplots(figsize=(6, 5))
    for cat, color in palette.items():
        mask = (meta[color_col] == cat).values
        if mask.sum() == 0:
            continue
        ax.scatter(coords[mask, 0], coords[mask, 1], s=30, c=color,
                   label=f"{cat} (n={int(mask.sum())})", alpha=0.85,
                   edgecolors="k", linewidths=0.3)
    ax.set_xlabel(f"PC1 ({var[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1]:.1f}%)")
    ax.set_title(title)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def boxplot(expr, gse):
    fig, ax = plt.subplots(figsize=(14, 4.5))
    data = [expr[c].values for c in expr.columns]
    ax.boxplot(data, showfliers=False, widths=0.6, patch_artist=True,
               boxprops=dict(facecolor="#cfe2f3", color="#3b6fa0"),
               medianprops=dict(color="#c0392b"))
    ax.set_xticks(range(1, len(expr.columns) + 1))
    ax.set_xticklabels(expr.columns, rotation=90)
    ax.set_ylabel("expression (log2)")
    ax.set_title(f"{gse} per-sample boxplot (n={expr.shape[1]})")
    fig.tight_layout()
    out = FIG / f"fig2d_boxplot_{gse.lower()}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def density(expr, gse):
    fig, ax = plt.subplots(figsize=(8, 5))
    for c in expr.columns:
        expr[c].plot(kind="density", ax=ax, color="grey", alpha=0.25, linewidth=0.6)
    ax.set_xlabel("expression (log2)")
    ax.set_ylabel("density")
    ax.set_title(f"{gse} per-sample density (n={expr.shape[1]})")
    fig.tight_layout()
    out = FIG / f"fig2e_density_{gse.lower()}.png"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def main():
    # GSE66360
    expr6, meta6 = load("GSE66360")
    pca_plot(expr6, meta6, "role", "GSE66360 PCA by role", FIG / "fig2a_pca_gse66360_by_role.png",
             {"train": "#1f77b4", "internal_val": "#ff7f0e"})
    pca_plot(expr6, meta6, "group", "GSE66360 PCA by case/control", FIG / "fig2b_pca_gse66360_by_group.png",
             {"case": "#d62728", "control": "#2ca02c"})
    boxplot(expr6, "GSE66360")
    density(expr6, "GSE66360")

    # GSE48060
    expr8, meta8 = load("GSE48060")
    pca_plot(expr8, meta8, "group", "GSE48060 PCA by case/control", FIG / "fig2c_pca_gse48060_by_group.png",
             {"case": "#d62728", "control": "#2ca02c"})

    print("[04_qc_plot] 完成")


if __name__ == "__main__":
    main()