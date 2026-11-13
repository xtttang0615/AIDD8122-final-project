# -*- coding: utf-8 -*-
"""16 P7 图表：ROC / PR / EF 曲线 / 消融条形图（300 dpi，字号>=8pt）。

输入：results/tables/ligand_oof_predictions.csv、ligand_model_perf.csv、ablation_summary.csv
输出：results/figures/fig4_ligand_roc.png、fig4b_ligand_pr.png、
      fig5_ef_curve.png、fig6_ablation.png
运行：C:/Anaconda/Anaconda3/python.exe src/16_ligand_figures.py
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
from sklearn.metrics import roc_curve, precision_recall_curve, roc_auc_score, average_precision_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import MAIN_TARGETS, ROOT, TAB  # noqa: E402

FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 100, "savefig.dpi": 300,
})
COLORS = {"PDE4B": "#1f77b4", "PFKFB3": "#2ca02c", "EPAS1": "#d62728"}


def plot_roc(preds, perf, main_model):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    for ax, target in zip(axes, MAIN_TARGETS):
        sub = preds[(preds["target"] == target) & (preds["model"] == main_model)]
        fpr, tpr, _ = roc_curve(sub["y"], sub["prob"])
        row = perf[(perf["target"] == target) & (perf["model"] == main_model)
                   & (perf["features"] == "ecfp4")].iloc[0]
        ax.plot(fpr, tpr, color=COLORS[target], lw=1.8,
                label=f"AUC={row['roc_auc']:.3f} ({row['ci_lo']:.3f}-{row['ci_hi']:.3f})")
        ax.plot([0, 1], [0, 1], ls="--", c="grey", lw=1)
        ax.set_title(f"{target} (n={len(sub)})")
        ax.set_xlabel("1 - specificity")
        ax.set_ylabel("sensitivity")
        ax.legend(loc="lower right", frameon=False)
    fig.suptitle(f"Ligand-based model ROC (main model: {main_model})")
    fig.tight_layout()
    fig.savefig(FIG / "fig4_ligand_roc.png", dpi=300)
    plt.close(fig)
    print("[fig] fig4_ligand_roc.png")


def plot_pr(preds, perf, main_model):
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 3.2))
    for ax, target in zip(axes, MAIN_TARGETS):
        sub = preds[(preds["target"] == target) & (preds["model"] == main_model)]
        prec, rec, _ = precision_recall_curve(sub["y"], sub["prob"])
        row = perf[(perf["target"] == target) & (perf["model"] == main_model)
                   & (perf["features"] == "ecfp4")].iloc[0]
        base = sub["y"].mean()
        ax.plot(rec, prec, color=COLORS[target], lw=1.8,
                label=f"AP={row['pr_auc']:.3f} ({row['pr_lo']:.3f}-{row['pr_hi']:.3f})")
        ax.axhline(base, ls="--", c="grey", lw=1, label=f"prevalence={base:.3f}")
        ax.set_title(f"{target} (n={len(sub)})", fontsize=10)
        ax.set_xlabel("recall", fontsize=9.5)
        ax.set_ylabel("precision", fontsize=9.5)
        ax.tick_params(labelsize=8.5)
        ax.legend(loc="lower left", bbox_to_anchor=(0.03, 0.12), frameon=False,
                  fontsize=7.5, borderpad=0.2, labelspacing=0.25, handlelength=1.6)
    fig.suptitle(f"Precision-recall curves (main model: {main_model})", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIG / "fig4b_ligand_pr.png", dpi=300)
    plt.close(fig)
    print("[fig] fig4b_ligand_pr.png")


def plot_ef(preds, main_model):
    fig, ax = plt.subplots(figsize=(6, 4.6))
    xs = np.linspace(0.005, 0.20, 40)
    for target in MAIN_TARGETS:
        sub = preds[(preds["target"] == target) & (preds["model"] == main_model)]
        y = sub["y"].to_numpy()
        p = sub["prob"].to_numpy()
        order = np.argsort(-p)
        base = y.mean()
        efs = []
        for x in xs:
            k = max(1, int(round(x * len(y))))
            efs.append((y[order[:k]].sum() / k) / base)
        ax.plot(xs * 100, efs, color=COLORS[target], lw=1.6, label=target)
    ax.axhline(1.0, ls="--", c="grey", lw=1, label="random")
    ax.set_xlabel("top % of ranked library")
    ax.set_ylabel("enrichment factor (EF)")
    ax.set_title(f"EF curve (main model: {main_model})")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "fig5_ef_curve.png", dpi=300)
    plt.close(fig)
    print("[fig] fig5_ef_curve.png")


def plot_ablation(ab):
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    specs = [("A1_features", "features", "A1: feature set"),
             ("A2_model", "model", "A2: model"),
             ("A3_threshold", "setting", "A3: pChEMBL threshold (PDE4B)")]
    for ax, (abl, col, title) in zip(axes, specs):
        sub = ab[ab["ablation"] == abl]
        if abl == "A3_threshold":
            r = sub.iloc[0]
            if pd.isna(r["pr_auc"]):
                # no-op：样本集合未变化，不绘制误导性比较
                ax.text(0.5, 0.5, "A3: no-op\n"
                        "pChEMBL>=7 did not change the PDE4B\npositive set (already capped at 800)",
                        ha="center", va="center", fontsize=8, wrap=True)
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_title(f"{title} (not applicable)")
                continue
            base_rf = ab[(ab["ablation"] == "A2_model") & (ab["target"] == "PDE4B")
                         & (ab["model"] == "rf")]
            base_val = base_rf["pr_auc"].iloc[0] if len(base_rf) else np.nan
            labels = [f">=6 (base)\n{base_val:.3f}", f">=7\n{sub['pr_auc'].iloc[0]:.3f}"]
            vals = [base_val, sub["pr_auc"].iloc[0]]
            ax.bar(labels, vals, color=["#7f7f7f", "#1f77b4"])
        else:
            piv = sub.pivot_table(index=col, values="pr_auc", aggfunc="mean")
            ax.bar(piv.index.astype(str), piv["pr_auc"].values,
                   color=plt.cm.tab10(np.arange(len(piv))))
            ax.set_ylabel("mean PR-AUC")
        ax.set_title(title)
        ax.tick_params(axis="x", rotation=20)
    fig.suptitle("Ablation summary (PR-AUC)")
    fig.tight_layout()
    fig.savefig(FIG / "fig6_ablation.png", dpi=300)
    plt.close(fig)
    print("[fig] fig6_ablation.png")


def main():
    preds = pd.read_csv(TAB / "ligand_oof_predictions.csv", encoding="utf-8-sig")
    perf = pd.read_csv(TAB / "ligand_model_perf.csv", encoding="utf-8-sig")
    ab = pd.read_csv(TAB / "ablation_summary.csv", encoding="utf-8-sig")
    main_model = perf.loc[perf["is_main"], "model"].iloc[0]
    plot_roc(preds, perf, main_model)
    plot_pr(preds, perf, main_model)
    plot_ef(preds, main_model)
    plot_ablation(ab)
    print("[16_ligand_figures] 完成")


if __name__ == "__main__":
    main()