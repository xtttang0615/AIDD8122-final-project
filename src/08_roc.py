# -*- coding: utf-8 -*-
"""P4 图：三层 ROC 与校准曲线（300 dpi，字号>=8pt）。

输入：results/tables/predictions.csv、results/tables/model_performance.csv
输出：results/figures/fig3_roc.png、results/figures/fig3b_calibration.png
运行：C:/Anaconda/Anaconda3/python.exe src/08_roc.py
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
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.calibration import calibration_curve

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
FIG = ROOT / "results" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 8, "ytick.labelsize": 8, "legend.fontsize": 8,
    "figure.dpi": 100, "savefig.dpi": 300,
})

COLORS = {"train_nested_full": "#1f77b4", "train_nested_fixedpool": "#17becf",
          "internal_val": "#2ca02c", "external": "#d62728", "external_zscored": "#ff7f0e"}
LABELS = {"train_nested_full": "train (fully nested CV)",
          "train_nested_fixedpool": "train (nested, fixed pool)",
          "internal_val": "internal_val (n=56)",
          "external": "external (GSE48060, n=52)",
          "external_zscored": "external z-scored (sensitivity)"}

CURVES = ["train_nested_full", "internal_val", "external", "external_zscored"]


def plot_roc(preds, perf):
    fig, ax = plt.subplots(figsize=(6, 5.5))
    for set_name in CURVES:
        sub = preds[preds["set"] == set_name]
        fpr, tpr, _ = roc_curve(sub["y_true"], sub["y_prob"])
        row = perf[perf["set"] == set_name].iloc[0]
        ax.plot(fpr, tpr, color=COLORS[set_name], lw=1.8,
                label=f"{LABELS[set_name]} AUC={row['auc']:.3f} "
                      f"({row['auc_lo']:.3f}-{row['auc_hi']:.3f})")
    ax.plot([0, 1], [0, 1], ls="--", color="grey", lw=1, label="chance")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.01)
    ax.set_xlabel("1 - specificity")
    ax.set_ylabel("sensitivity")
    ax.set_title("AMI diagnostic signature: ROC (three-tier)")
    ax.legend(loc="lower right", frameon=False)
    out = FIG / "fig3_roc.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def plot_calibration(preds):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    ax.plot([0, 1], [0, 1], ls="--", color="grey", lw=1, label="perfect")
    for set_name in ["external", "internal_val"]:
        sub = preds[preds["set"] == set_name]
        n_bins = 5 if len(sub) < 80 else 8
        frac_pos, mean_pred = calibration_curve(sub["y_true"], sub["y_prob"],
                                                n_bins=n_bins, strategy="quantile")
        brier = float(np.mean((sub["y_prob"] - sub["y_true"]) ** 2))
        ax.plot(mean_pred, frac_pos, marker="o", ms=4, lw=1.5, color=COLORS[set_name],
                label=f"{LABELS[set_name]}  Brier={brier:.3f}")
    ax.set_xlabel("mean predicted probability")
    ax.set_ylabel("observed fraction of cases")
    ax.set_title("Calibration (external primary)")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc="lower right", frameon=False)
    out = FIG / "fig3b_calibration.png"
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"[fig] {out.name}")


def main():
    preds = pd.read_csv(TAB / "predictions.csv", encoding="utf-8-sig")
    perf = pd.read_csv(TAB / "model_performance.csv", encoding="utf-8-sig")
    plot_roc(preds, perf)
    plot_calibration(preds)

    with open(ROOT / "logs" / "model_log.md", "a", encoding="utf-8") as f:
        f.write("\n## 图（08）\n")
        f.write("- fig3_roc.png：三条 ROC 同图，标 AUC 与 95% CI，对角虚线为 chance\n")
        f.write("- fig3b_calibration.png：校准曲线，external 为主，另附 internal_val 与 Brier 分数\n")
    print("[08_roc] 完成")


if __name__ == "__main__":
    main()