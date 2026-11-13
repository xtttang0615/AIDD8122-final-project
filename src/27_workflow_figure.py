# -*- coding: utf-8 -*-
"""27 重绘工作流程图（v2）：修正列间距/箭头，补回被截断的 nested CV，
并按成品尺寸设计字号（图宽 6.5 in，插入 Word 后约 6.3 in，字号≈7 pt）。

输出：workflow_pipeline.png（项目根，400 dpi）
运行：C:/Anaconda/Anaconda3/python.exe src/27_workflow_figure.py
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "savefig.dpi": 400})

W, H = 6.5, 3.95           # inches — final print size
COLW = 0.208
XS = [0.012, 0.268, 0.524, 0.780]
ROWY = [0.775, 0.615, 0.455, 0.295]
BOXH = 0.115
TITLE_Y = 0.935

GROUPS = [
    ("Transcriptomic branch", "#dbe9f6", [
        "GEO download\nGSE66360 / GSE48060\n+ GPL570 annotation",
        "Parse + QC\n(log2 check, probe→gene,\nmax-mean probe)",
        "limma DE\n(training samples only;\n|logFC|>1 & FDR<0.05)",
        "Signature\nbootstrap 200×\n(RF top20 ∩ LASSO)",
    ]),
    ("Target annotation", "#e8f4d9", [
        "Open Targets\nMI association",
        "STRING\ninteraction degree",
        "ChEMBL\nligands / max pChEMBL",
        "Rule-based screen\n(strict / relaxed)",
    ]),
    ("Ligand benchmark", "#fde8d7", [
        "Actives + decoys\n(MolStandardize,\nInChIKey14 dedup)",
        "ECFP4 features\n+ 9 descriptors",
        "Scaffold splits\n(StratifiedGroupKFold)",
        "RF / LR / XGB\nOOF evaluation\n+ hard negatives",
    ]),
    ("Independent evaluation", "#f3e0f0", [
        "Approved drug pool\n(max_phase = 4,\nexclude train overlap)",
        "Score 3 targets\n→ consensus",
        "MI-indication set\n(exploratory only)",
        "Shortlist\n+ model-score ranking",
    ]),
]


def box(ax, x, y, w, h, text, fc, fs=7.2):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                                boxstyle="round,pad=0.004,rounding_size=0.012",
                                linewidth=0.8, edgecolor="#4a4a4a", facecolor=fc))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fs, linespacing=1.25)


def arrow(ax, p1, p2, dashed=False):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle="-|>", mutation_scale=7,
                                 linewidth=0.9, color="#555555",
                                 linestyle="--" if dashed else "-",
                                 shrinkA=0, shrinkB=0))


def main():
    fig, ax = plt.subplots(figsize=(W, H))
    ax.set_xlim(-0.004, 1.004)
    ax.set_ylim(0, 1)
    ax.axis("off")

    for gi, (title, color, items) in enumerate(GROUPS):
        x0 = XS[gi]
        ax.text(x0 + COLW / 2, TITLE_Y, title, ha="center", va="center",
                fontsize=7.1, weight="bold", color="#222222")
        for i, text in enumerate(items):
            box(ax, x0, ROWY[i], COLW, BOXH, text, color)
            if i:
                arrow(ax, (x0 + COLW / 2, ROWY[i - 1]),
                      (x0 + COLW / 2, ROWY[i] + BOXH))

    # inter-column flow, at the level of the last step of each column
    ymid = ROWY[3] + BOXH / 2
    for gi in range(3):
        arrow(ax, (XS[gi] + COLW, ymid), (XS[gi + 1], ymid))

    # nested cross-validation: spans branch 1-2, validates the signature
    yb = 0.125
    ncv_w = XS[1] + COLW - XS[0]
    box(ax, XS[0], yb, ncv_w, 0.100,
        "Nested CV — fold-wise limma + selection + L2 LR,\n"
        "the whole chain re-run inside each training fold (AUC 0.961)",
        "#eeeeee", fs=6.9)
    arrow(ax, (XS[0] + 0.055, yb + 0.100), (XS[0] + 0.055, ROWY[3]), dashed=True)

    ax.text(0.5, 0.035,
            "Executed scope only: no GNN, no DeepPurpose, no molecular docking or free-energy calculation.\n"
            "Model outputs are uncalibrated model scores, not binding probabilities.",
            ha="center", va="center", fontsize=5.9, style="italic", color="#555555")

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    fig.savefig(ROOT / "workflow_pipeline.png", dpi=400, facecolor="white")
    plt.close(fig)
    print("[27] workflow_pipeline.png 已更新 (v2)")


if __name__ == "__main__":
    main()
