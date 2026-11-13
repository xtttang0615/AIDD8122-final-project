# -*- coding: utf-8 -*-
"""DEG 图：火山图、前 50 热图、top2000 PCA、离群检查（300 dpi，字号>=8pt）。

输入：results/tables/deg_train.csv、data/processed/GSE66360_gene.csv、data/processed/GSE66360_meta.csv
输出：results/figures/fig2g_volcano.png、fig2h_heatmap_top50.png、fig2f_pca_top2000.png
       logs/deg_log.md（追加）
运行：C:/Anaconda/Anaconda3/python.exe src/06_deg_plots.py
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
from adjustText import adjust_text

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
LOGS = ROOT / "logs"
PROCESSED = ROOT / "data" / "processed"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 100, "savefig.dpi": 300,
})


def load():
    expr = pd.read_csv(PROCESSED / "GSE66360_gene.csv", index_col=0, encoding="utf-8-sig")
    meta = pd.read_csv(PROCESSED / "GSE66360_meta.csv", encoding="utf-8-sig").set_index("sample_id")
    meta = meta.loc[expr.columns]
    return expr, meta


def volcano(deg):
    up = (deg["logFC"] > 1) & (deg["adj.P.Val"] < 0.05)
    down = (deg["logFC"] < -1) & (deg["adj.P.Val"] < 0.05)
    ns = ~(up | down)
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(deg.loc[ns, "logFC"], -np.log10(deg.loc[ns, "adj.P.Val"]),
               s=5, c="lightgrey", alpha=0.5, label="ns")
    ax.scatter(deg.loc[up, "logFC"], -np.log10(deg.loc[up, "adj.P.Val"]),
               s=10, c="#d62728", alpha=0.85, label=f"up ({int(up.sum())})")
    ax.scatter(deg.loc[down, "logFC"], -np.log10(deg.loc[down, "adj.P.Val"]),
               s=10, c="#1f77b4", alpha=0.85, label=f"down ({int(down.sum())})")
    ax.axhline(-np.log10(0.05), color="k", ls="--", lw=1)
    ax.axvline(1, color="k", ls="--", lw=1)
    ax.axvline(-1, color="k", ls="--", lw=1)
    ax.set_xlabel("logFC")
    ax.set_ylabel("-log10(adj.P.Val)")
    ax.set_title("GSE66360 train volcano")
    ax.legend(frameon=False)

    up_t = deg[up].nlargest(10, "logFC")
    dn_t = deg[down].nsmallest(10, "logFC")
    combo = pd.concat([up_t, dn_t])
    texts = []
    for _, r in combo.iterrows():
        texts.append(ax.text(r["logFC"], -np.log10(r["adj.P.Val"]), r["gene"], fontsize=8))
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5),
                only_move={"points": "y", "texts": "y"},
                expand_points=(1.2, 1.2), expand_text=(1.2, 1.2))
    ylim = ax.get_ylim()
    for t in texts:
        x, y = t.get_position()
        t.set_position((x, min(max(y, ylim[0]), ylim[1])))

    out = FIG / "fig2g_volcano.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def heatmap(expr, meta, deg):
    sig = deg[deg["sig"]]
    top50 = sig.nsmallest(50, "adj.P.Val").sort_values("logFC", ascending=False)
    genes = top50["gene"].tolist()
    mat = expr.loc[genes]
    train_ids = [s for s in meta.index if meta.loc[s, "role"] == "train" and s in mat.columns]
    groups = meta.loc[train_ids, "group"]
    order = [s for s in train_ids if groups[s] == "control"] + [s for s in train_ids if groups[s] == "case"]
    sub = mat[order]
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, np.nan), axis=0).fillna(0)

    fig, ax = plt.subplots(figsize=(12, 8))
    im = ax.imshow(z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=8)
    ax.set_title("GSE66360 train top50 DEG heatmap (row z-score)")
    cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("z-score")

    out = FIG / "fig2h_heatmap_top50.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def pca_top2000(expr, meta):
    var = expr.var(axis=1)
    top = var.nlargest(2000).index
    X = expr.loc[top].T.values
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)
    var_r = pca.explained_variance_ratio_ * 100

    fig, ax = plt.subplots(figsize=(6, 5))
    for cat, col in {"case": "#d62728", "control": "#2ca02c"}.items():
        mask = (meta["group"] == cat).values
        ax.scatter(coords[mask, 0], coords[mask, 1], s=30, c=col, alpha=0.85,
                   edgecolors="k", linewidths=0.3, label=f"{cat} (n={int(mask.sum())})")
    ax.set_xlabel(f"PC1 ({var_r[0]:.1f}%)")
    ax.set_ylabel(f"PC2 ({var_r[1]:.1f}%)")
    ax.set_title("GSE66360 PCA (top-2000 variable genes)")
    ax.legend(frameon=False)

    out = FIG / "fig2f_pca_top2000.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")

    z1 = (coords[:, 0] - coords[:, 0].mean()) / coords[:, 0].std()
    outliers = [(sid, round(float(zv), 2)) for sid, zv in zip(expr.columns, z1) if abs(zv) > 3]
    return outliers


def main():
    expr, meta = load()
    deg = pd.read_csv(TAB / "deg_train.csv", encoding="utf-8-sig")
    deg["sig"] = deg["sig"].astype(str).str.upper().isin(["TRUE", "1", "T"])

    volcano(deg)
    heatmap(expr, meta, deg)
    outliers = pca_top2000(expr, meta)

    with open(LOGS / "deg_log.md", "a", encoding="utf-8") as f:
        f.write("\n## 图与离群检查（06）\n")
        f.write("- 火山图：fig2g_volcano.png\n")
        f.write("- 热图前50：fig2h_heatmap_top50.png\n")
        f.write("- PCA top2000：fig2f_pca_top2000.png\n")
        f.write(f"- PC1 |z|>3 离群样本：{outliers if outliers else '无'}\n")

    print(f"[06_deg_plots] 离群样本(|z|>3)：{outliers if outliers else '无'}")


if __name__ == "__main__":
    main()
