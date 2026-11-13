# -*- coding: utf-8 -*-
"""12 配体基打分模型：LR / RF / XGB，骨架切分 OOF 评估。

输入：data/processed/ligand_features_ecfp4.npz、ligand_features_desc.csv、scaffold_splits.csv
输出：results/tables/ligand_model_perf.csv、ligand_model_byfold.csv、ligand_oof_predictions.csv
      logs/p7_experiments.md
运行：C:/Anaconda/Anaconda3/python.exe src/12_ligand_model.py

主模型预注册规则（先写日志再跑）：默认 XGBoost；若三靶点 mean PR-AUC 上 XGB 与 RF 差 < 0.02，
取 RF。选定后不再更换。切分来自 11_splits 固化表，本脚本不重新生成。
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
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import (DESC_COLS, LOGS, MAIN_TARGETS, SEED, TAB,  # noqa: E402
                          ALL_TARGETS, bootstrap_ci, load_features, load_splits,
                          make_model, metrics)

MODELS = ["lr", "rf", "xgb"]


def build_matrix(ids, ecfp, desc, mode):
    Xf = np.vstack([ecfp[m] for m in ids]).astype(np.float32)
    if mode == "ecfp4":
        return Xf
    D = desc.loc[ids, DESC_COLS].to_numpy(dtype=np.float32)
    if mode == "desc":
        return D
    return np.hstack([Xf, D])


def run_cv(df, ecfp, desc, model_name, mode):
    """按固化 fold 跑 OOF，描述符在每个训练折内 z-score。"""
    ids = df["molecule_chembl_id"].tolist()
    y = df["y"].to_numpy()
    folds = df["fold"].to_numpy()
    X = build_matrix(ids, ecfp, desc, mode)
    oof = np.full(len(df), np.nan)
    by_fold = []

    holdout = df["scheme"].iloc[0] == "holdout80_20"
    # 单次 hold-out 只跑一次，避免 fold=0/1 各拟合一遍造成重复记录
    fold_iter = [1] if holdout else sorted(np.unique(folds))

    for k in fold_iter:
        if holdout:
            tr, te = folds == 0, folds == 1
        else:
            tr, te = folds != k, folds == k
        if te.sum() == 0 or len(np.unique(y[tr])) < 2:
            continue

        Xtr, Xte = X[tr], X[te]
        if mode in ("desc", "both"):
            sc = StandardScaler()
            if mode == "desc":
                Xtr = sc.fit_transform(Xtr)
                Xte = sc.transform(Xte)
            else:
                Xtr = Xtr.copy()
                Xte = Xte.copy()
                Xtr[:, -len(DESC_COLS):] = sc.fit_transform(Xtr[:, -len(DESC_COLS):])
                Xte[:, -len(DESC_COLS):] = sc.transform(Xte[:, -len(DESC_COLS):])

        model = make_model(model_name, y[tr])
        model.fit(Xtr, y[tr])
        p = model.predict_proba(Xte)[:, 1]
        oof[te] = p

        if len(np.unique(y[te])) > 1:
            from sklearn.metrics import roc_auc_score, average_precision_score
            by_fold.append({"fold": int(k), "n_pos": int(y[te].sum()),
                            "n_neg": int((1 - y[te]).sum()),
                            "pr_auc": float(average_precision_score(y[te], p)),
                            "roc_auc": float(roc_auc_score(y[te], p))})
        else:
            by_fold.append({"fold": int(k), "n_pos": int(y[te].sum()),
                            "n_neg": int((1 - y[te]).sum()),
                            "pr_auc": np.nan, "roc_auc": np.nan})
    return oof, by_fold


def main():
    ecfp, desc = load_features()
    splits = load_splits()

    exp_log = ["# P7 实验记录", "",
               f"## {datetime.datetime.now().isoformat(timespec='seconds')} 预注册主模型规则",
               "- 默认主模型：XGBoost",
               "- 规则：若三主靶点（PDE4B/PFKFB3/EPAS1）mean PR-AUC 上 XGB 与 RF 之差 < 0.02，取 RF（更简单者胜）",
               "- 该规则在看结果之前写下，选定后不再更换",
               ""]

    perf_rows, byfold_rows, pred_rows = [], [], []
    summary = {}

    for target in ALL_TARGETS:
        df = splits[splits["target"] == target].reset_index(drop=True)
        y = df["y"].to_numpy()
        tag = "（exploratory）" if target == "PDK4" else ""
        print(f"\n=== {target} n={len(df)} pos={int(y.sum())} {tag}")

        # 随机基线
        rng = np.random.default_rng(SEED)
        p_rand = rng.random(len(df))
        m_rand = metrics(y, p_rand)
        perf_rows.append({"target": target, "model": "random", "features": "-",
                          "roc_auc": m_rand["roc_auc"], "ci_lo": np.nan, "ci_hi": np.nan,
                          "pr_auc": m_rand["pr_auc"], "pr_lo": np.nan, "pr_hi": np.nan,
                          "ef1": m_rand["ef1"], "ef5": m_rand["ef5"], "p@10": m_rand["p@10"],
                          "is_main": False, "exploratory": target == "PDK4"})

        for model_name in MODELS:
            oof, bf = run_cv(df, ecfp, desc, model_name, "ecfp4")
            mask = ~np.isnan(oof)          # PDK4 只评估留出集
            m = metrics(y[mask], oof[mask])
            lo, hi, plo, phi = bootstrap_ci(y[mask], oof[mask])
            perf_rows.append({"target": target, "model": model_name, "features": "ecfp4",
                              "roc_auc": m["roc_auc"], "ci_lo": lo, "ci_hi": hi,
                              "pr_auc": m["pr_auc"], "pr_lo": plo, "pr_hi": phi,
                              "ef1": m["ef1"], "ef5": m["ef5"], "p@10": m["p@10"],
                              "is_main": False, "exploratory": target == "PDK4"})
            for b in bf:
                byfold_rows.append({"target": target, "model": model_name, **b})
            for mid, yy, pp in zip(df["molecule_chembl_id"][mask], y[mask], oof[mask]):
                pred_rows.append({"target": target, "model": model_name,
                                  "molecule_chembl_id": mid, "y": int(yy), "prob": float(pp)})
            summary[(target, model_name)] = m["pr_auc"]
            print(f"  {model_name}: ROC={m['roc_auc']:.4f} PR={m['pr_auc']:.4f} "
                  f"EF1={m['ef1']:.2f} EF5={m['ef5']:.2f} p@10={m['p@10']:.3f}")

        # 仅描述符基线（用 LR，代表线性可用的物化信息）
        oof_d, _ = run_cv(df, ecfp, desc, "lr", "desc")
        mask_d = ~np.isnan(oof_d)
        m_d = metrics(y[mask_d], oof_d[mask_d])
        perf_rows.append({"target": target, "model": "lr", "features": "desc",
                          "roc_auc": m_d["roc_auc"], "ci_lo": np.nan, "ci_hi": np.nan,
                          "pr_auc": m_d["pr_auc"], "pr_lo": np.nan, "pr_hi": np.nan,
                          "ef1": m_d["ef1"], "ef5": m_d["ef5"], "p@10": m_d["p@10"],
                          "is_main": False, "exploratory": target == "PDK4"})
        print(f"  lr/desc 基线：ROC={m_d['roc_auc']:.4f} PR={m_d['pr_auc']:.4f}")

    perf = pd.DataFrame(perf_rows)
    byfold = pd.DataFrame(byfold_rows)

    # ---- 主模型选择（预注册规则）
    mean_pr = {}
    for mn in ("xgb", "rf"):
        vals = [summary[(t, mn)] for t in MAIN_TARGETS]
        mean_pr[mn] = float(np.mean(vals))
    diff = abs(mean_pr["xgb"] - mean_pr["rf"])
    main_model = "rf" if diff < 0.02 else "xgb"
    exp_log += [f"## {datetime.datetime.now().isoformat(timespec='seconds')} 主模型选定",
                f"- 三靶点 mean PR-AUC：XGB={mean_pr['xgb']:.4f}，RF={mean_pr['rf']:.4f}，"
                f"差值={diff:.4f}（阈值 0.02）",
                f"- 结果：主模型 = **{main_model}**",
                ""]

    perf.loc[(perf["model"] == main_model) & (perf["features"] == "ecfp4"), "is_main"] = True
    perf = perf.sort_values(["target", "is_main", "pr_auc"], ascending=[True, False, False])
    perf.to_csv(TAB / "ligand_model_perf.csv", index=False, encoding="utf-8-sig")
    byfold.to_csv(TAB / "ligand_model_byfold.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(pred_rows).to_csv(TAB / "ligand_oof_predictions.csv", index=False,
                                   encoding="utf-8-sig")

    # 折间 PR-AUC 波动
    exp_log.append("## 折间 PR-AUC（mean±SD，仅稳定性参考，不作为判据）")
    for target in ALL_TARGETS:
        sub = byfold[(byfold["target"] == target) & (byfold["model"] == main_model)]
        vals = sub["pr_auc"].dropna()
        exp_log.append(f"- {target}（{main_model}）：{vals.mean():.4f} ± {vals.std(ddof=0):.4f}"
                       f"（n_folds={len(vals)}）")
    exp_log += ["", "## 运行记录（时间/配置/结果/是否主结果）"]
    for r in perf.itertuples(index=False):
        exp_log.append(f"- {datetime.datetime.now().isoformat(timespec='seconds')} | "
                       f"{r.target} | {r.model} | {r.features} | "
                       f"ROC={r.roc_auc:.4f} PR={r.pr_auc:.4f} | "
                       f"{'主结果' if r.is_main else '对照'}")
    exp_log.append("")
    (LOGS / "p7_experiments.md").write_text("\n".join(exp_log), encoding="utf-8")

    print(f"\n主模型：{main_model}（XGB mean PR={mean_pr['xgb']:.4f} vs RF={mean_pr['rf']:.4f}）")
    print("[12_ligand_model] 完成")


if __name__ == "__main__":
    main()