# -*- coding: utf-8 -*-
"""13 消融实验：特征 / 模型 / 活性阈值敏感性（复用同一份固化切分）。

输入：data/processed/ligand_features_ecfp4.npz、ligand_features_desc.csv、scaffold_splits.csv
      results/tables/ligand_model_perf.csv（读取主模型名）
输出：results/tables/ablation_summary.csv、logs/p7_experiments.md（追加）
运行：C:/Anaconda/Anaconda3/python.exe src/13_ablation.py

A1 特征：ECFP4 / 仅描述符 / ECFP4+描述符（主模型）
A2 模型：LR / RF / XGB（ECFP4）
A3 阈值：PDE4B 正样本 pChEMBL>=7 重建，诱饵不变，同一切分下对比 PR-AUC
"""
from pathlib import Path
import sys
import datetime

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import (LOGS, MAIN_TARGETS, TAB, load_features, load_splits,  # noqa: E402
                          metrics)

# 直接复用 12 号脚本中的 run_cv：以文件方式加载
import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "ligand_model", Path(__file__).resolve().parent / "12_ligand_model.py")
_lm = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_lm)

A3_THRESH = 7.0


def main():
    ecfp, desc = load_features()
    splits = load_splits()
    perf = pd.read_csv(TAB / "ligand_model_perf.csv", encoding="utf-8-sig")
    main_model = perf.loc[perf["is_main"], "model"].iloc[0]
    print(f"主模型（来自 12 号）：{main_model}")

    rows = []

    # ---- A1 特征
    for feat in ["ecfp4", "desc", "both"]:
        for target in MAIN_TARGETS:
            df = splits[splits["target"] == target].reset_index(drop=True)
            oof, _ = _lm.run_cv(df, ecfp, desc, main_model, feat)
            m = metrics(df["y"].to_numpy(), oof)
            rows.append({"ablation": "A1_features", "target": target, "model": main_model,
                         "features": feat, "setting": "default", "pr_auc": m["pr_auc"],
                         "roc_auc": m["roc_auc"], "ef1": m["ef1"], "note": ""})
            print(f"A1 {target} {feat}: PR={m['pr_auc']:.4f}")

    # ---- A2 模型
    for mn in ["lr", "rf", "xgb"]:
        for target in MAIN_TARGETS:
            df = splits[splits["target"] == target].reset_index(drop=True)
            oof, _ = _lm.run_cv(df, ecfp, desc, mn, "ecfp4")
            m = metrics(df["y"].to_numpy(), oof)
            rows.append({"ablation": "A2_model", "target": target, "model": mn,
                         "features": "ecfp4", "setting": "default", "pr_auc": m["pr_auc"],
                         "roc_auc": m["roc_auc"], "ef1": m["ef1"], "note": ""})
            print(f"A2 {target} {mn}: PR={m['pr_auc']:.4f}")

    # ---- A3 阈值敏感性（PDE4B）
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    pde = splits[splits["target"] == "PDE4B"].copy()
    pos_ids = set(act[(act["target"] == "PDE4B") & (act["pchembl"] >= A3_THRESH)]
                  ["molecule_chembl_id"])
    a3 = pde[(pde["y"] == 0) | (pde["molecule_chembl_id"].isin(pos_ids))].reset_index(drop=True)
    n_pos, n_neg = int(a3["y"].sum()), int((1 - a3["y"]).sum())
    n_pos_base = int((splits[(splits.target == "PDE4B") & (splits.y == 1)]).shape[0])
    base = perf[(perf["target"] == "PDE4B") & (perf["is_main"])].iloc[0]

    if n_pos == n_pos_base:
        # 样本集合未变化，不构成有效敏感性分析
        rows.append({"ablation": "A3_threshold", "target": "PDE4B", "model": main_model,
                     "features": "ecfp4", "setting": f"pchembl>={A3_THRESH:.0f}",
                     "pr_auc": np.nan, "roc_auc": np.nan, "ef1": np.nan,
                     "note": "no-op：正样本数未变化（PDE4B 已按 pChEMBL 截断取前 800，"
                             "最小值 7.58），该阈值不构成有效敏感性分析，不作比较"})
        print(f"A3 PDE4B pChEMBL>=7: 正样本 {n_pos} 与基线 {n_pos_base} 相同 -> no-op，已标记")
    else:
        oof, _ = _lm.run_cv(a3, ecfp, desc, main_model, "ecfp4")
        m = metrics(a3["y"].to_numpy(), oof)
        note = (f"仅比较同一切分下的相对变化，正负比 800:6238 -> {n_pos}:{n_neg}，"
                "正负比已改变，不作绝对性能解读")
        rows.append({"ablation": "A3_threshold", "target": "PDE4B", "model": main_model,
                     "features": "ecfp4", "setting": f"pchembl>={A3_THRESH:.0f}",
                     "pr_auc": m["pr_auc"], "roc_auc": m["roc_auc"], "ef1": m["ef1"], "note": note})
        print(f"A3 PDE4B pChEMBL>=7: n_pos={n_pos} n_neg={n_neg} PR={m['pr_auc']:.4f} "
              f"(baseline PR={base['pr_auc']:.4f})")

    ab = pd.DataFrame(rows)
    ab.to_csv(TAB / "ablation_summary.csv", index=False, encoding="utf-8-sig")

    with open(LOGS / "p7_experiments.md", "a", encoding="utf-8") as f:
        f.write(f"\n## {datetime.datetime.now().isoformat(timespec='seconds')} 消融（13）\n")
        f.write(f"- 主模型：{main_model}；切分复用 scaffold_splits.csv，未重新生成\n")
        for r in ab.itertuples(index=False):
            f.write(f"- {r.ablation} | {r.target} | {r.model} | {r.features} | {r.setting} | "
                    f"PR={r.pr_auc:.4f} ROC={r.roc_auc:.4f} EF1={r.ef1:.2f}\n")
        f.write("- A3 note：仅比较同一切分下的相对变化，正负比已改变，不作绝对性能解读\n\n")

    print("[13_ablation] 完成")


if __name__ == "__main__":
    main()