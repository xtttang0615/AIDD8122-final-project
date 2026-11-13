# -*- coding: utf-8 -*-
"""P7d 基准可比性修正：加 prevalence/lift，PDE4B easy 下采样对齐。

输入：results/tables/ligand_model_hard_benchmark.csv、ligand_oof_predictions.csv
输出：results/tables/ligand_model_hard_benchmark.csv（新增列与 easy_downsampled 行）、
      logs/p7c_hard_benchmark.md（追加修正说明）
运行：C:/Anaconda/Anaconda3/python.exe src/19_benchmark_normalize.py

只重算，不重训模型、不调参、不改 ligand_model_perf.csv 与 scaffold_splits.csv。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, SEED, TAB, ef_at, precision_at  # noqa: E402

N_REPEAT = 20
TARGET_NEG = {"PDE4B": 1034}


def main():
    bench = pd.read_csv(TAB / "ligand_model_hard_benchmark.csv", encoding="utf-8-sig")

    # ---- 1. prevalence / lift
    bench["prevalence"] = bench["n_pos"] / (bench["n_pos"] + bench["n_neg"])
    bench["lift"] = bench["pr_auc"] / bench["prevalence"]
    bench["pr_sd"] = np.nan

    # ---- 2. PDE4B easy 下采样对齐（复用 12 号 OOF，不重训）
    oof = pd.read_csv(TAB / "ligand_oof_predictions.csv", encoding="utf-8-sig")
    sub = oof[(oof["target"] == "PDE4B") & (oof["model"] == "rf")]
    y = sub["y"].to_numpy()
    p = sub["prob"].to_numpy()
    pos_idx = np.where(y == 1)[0]
    neg_idx = np.where(y == 0)[0]
    n_neg_use = TARGET_NEG["PDE4B"]

    rng = np.random.default_rng(SEED)
    prs, rocs, ef1s, p10s = [], [], [], []
    for _ in range(N_REPEAT):
        pick = rng.choice(neg_idx, size=n_neg_use, replace=False)
        idx = np.concatenate([pos_idx, pick])
        yy, pp = y[idx], p[idx]
        prs.append(average_precision_score(yy, pp))
        rocs.append(roc_auc_score(yy, pp))
        ef1s.append(ef_at(yy, pp, 0.01))
        p10s.append(precision_at(yy, pp, 10))

    pr_mean, pr_sd = float(np.mean(prs)), float(np.std(prs, ddof=1))
    roc_mean = float(np.mean(rocs))
    prev = len(pos_idx) / (len(pos_idx) + n_neg_use)

    new_row = pd.DataFrame([{
        "target": "PDE4B", "benchmark": "easy_downsampled",
        "n_pos": int(len(pos_idx)), "n_neg": int(n_neg_use),
        "pr_auc": pr_mean, "roc_auc": roc_mean,
        "ef1": float(np.mean(ef1s)), "p@10": float(np.mean(p10s)),
        "prevalence": prev, "lift": pr_mean / prev, "pr_sd": pr_sd,
    }])
    bench = pd.concat([bench, new_row], ignore_index=True)
    bench = bench[["target", "benchmark", "n_pos", "n_neg", "prevalence", "pr_auc", "pr_sd",
                   "lift", "roc_auc", "ef1", "p@10"]]
    bench.to_csv(TAB / "ligand_model_hard_benchmark.csv", index=False, encoding="utf-8-sig")

    print(bench.round(4).to_string(index=False))
    print(f"PDE4B easy_downsampled：PR-AUC = {pr_mean:.4f} ± {pr_sd:.4f}（{N_REPEAT} 次）")

    # ---- 3/4/5 日志修正
    def get(t, b, col):
        r = bench[(bench.target == t) & (bench.benchmark == b)].iloc[0]
        return r[col]

    log_add = [
        "",
        "## P7d 可比性修正",
        "",
        "- 新增列：prevalence = n_pos/(n_pos+n_neg)；lift = PR-AUC / prevalence"
        "（随机排序的 PR-AUC 期望值即正例率）",
        f"- PDE4B easy 下采样对齐：从 6238 个负样本中随机抽 {n_neg_use} 个（seed={SEED}），"
        f"用已有 OOF 预测重算 {N_REPEAT} 次，PR-AUC = {pr_mean:.4f} ± {pr_sd:.4f}",
        "",
        "| target | benchmark | prevalence | PR-AUC | lift | ROC-AUC |",
        "|---|---|---|---|---|---|",
    ]
    for r in bench.itertuples(index=False):
        log_add.append(f"| {r.target} | {r.benchmark} | {r.prevalence:.4f} | "
                       f"{r.pr_auc:.4f} | {r.lift:.2f} | {r.roc_auc:.4f} |")

    log_add += [
        "",
        "ROC-AUC 对类别比例不敏感（0.9981 vs 0.9979），可作为跨基准比较的指标；"
        "PR-AUC 因正例率变化不可跨基准直接比较，须用 lift 或下采样对齐后比较。",
        "",
        "### 结论（统一表述）",
        "",
        "PFKFB3/EPAS1 正例率接近，hard 基准 PR-AUC 基本持平"
        f"（lift {get('PFKFB3','easy','lift'):.1f}→{get('PFKFB3','hard','lift'):.1f}、"
        f"{get('EPAS1','easy','lift'):.1f}→{get('EPAS1','hard','lift'):.1f}），"
        "提示存在靶点特异判别；PDE4B 因正例率由 11.4% 升至 43.6%，其 PR-AUC 不可直接比较，"
        f"下采样对齐后为 {pr_mean:.4f}（lift {get('PDE4B','easy','lift'):.1f}→"
        f"{get('PDE4B','easy_downsampled','lift'):.1f}）。",
        "",
        "三分层排序：本靶点配体 > 其他靶点配体 > 随机小分子；"
        "两层基准均接近饱和，仅支持相对排序用途。",
        "",
    ]
    with open(LOGS / "p7c_hard_benchmark.md", "a", encoding="utf-8") as f:
        f.write("\n".join(log_add))
    print("[19_benchmark_normalize] 完成")


if __name__ == "__main__":
    main()