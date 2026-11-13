# -*- coding: utf-8 -*-
"""18b 硬负样本基准（修订版）：先合并集合，再统一骨架五折。

输入：pool_actives.csv、pool_decoys.csv、scaffold_splits.csv、ligand_features_ecfp4.npz
输出：data/processed/hard_benchmark_splits.csv（新切分，原 scaffold_splits.csv 保留不动）
      results/tables/ligand_model_hard_benchmark.csv（覆盖前先归档为 _legacy）
      logs/p7c_hard_benchmark.md（追加修订记录）
运行：C:/Anaconda/Anaconda3/python.exe src/18b_hard_benchmark_v2.py

修订要点：原实现直接沿用各靶点独立生成的折号，合并后可能存在同一骨架跨折。
本版先构建完整困难基准，再统一按 Bemis-Murcko 骨架做 StratifiedGroupKFold(5, seed=42)，
并显式检查跨折冲突。模型与特征参数不变（RF500 + ECFP4）。
"""
from pathlib import Path
import shutil
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import (LOGS, MAIN_TARGETS, PROC, SEED, TAB, ef_at, metrics,  # noqa: E402
                          precision_at, load_features, load_splits)

T0 = time.time()


def main():
    ecfp, _ = load_features()
    splits = load_splits()
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    perf = pd.read_csv(TAB / "ligand_model_perf.csv", encoding="utf-8-sig")
    rng = np.random.default_rng(SEED)

    params = dict(n_estimators=500, class_weight="balanced_subsample",
                  random_state=SEED, n_jobs=-1)

    all_rows, split_rows, log = [], [], ["## P7c 修订（18b）：合并后统一骨架五折", ""]

    for target in MAIN_TARGETS:
        pos = splits[(splits["target"] == target) & (splits["y"] == 1)]
        pos_scaf = set(pos["scaffold"])
        pos_ik = set(act.loc[act["target"] == target, "inchikey14"])

        hard, dropped = [], 0
        for donor in [t for t in MAIN_TARGETS if t != target]:
            cand = splits[(splits["target"] == donor) & (splits["y"] == 1)]
            donor_ik = act[act["target"] == donor].set_index("molecule_chembl_id")["inchikey14"]
            for r in cand.itertuples(index=False):
                if donor_ik.get(r.molecule_chembl_id) in pos_ik:
                    dropped += 1
                    continue
                if r.scaffold in pos_scaf:
                    dropped += 1
                    continue
                hard.append({"molecule_chembl_id": r.molecule_chembl_id, "y": 0,
                             "scaffold": r.scaffold})

        easy = splits[(splits["target"] == target) & (splits["y"] == 0)]
        n_easy = min(len(easy), len(hard))
        easy_use = easy.iloc[rng.choice(len(easy), size=n_easy, replace=False)]

        df = pd.concat([
            pos[["molecule_chembl_id", "y", "scaffold"]],
            easy_use[["molecule_chembl_id", "y", "scaffold"]],
            pd.DataFrame(hard),
        ], ignore_index=True)

        # 标签冲突检查：同一 InChIKey14 同时为 0/1
        ik = pd.concat([act[["molecule_chembl_id", "inchikey14"]],
                        pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
                        [["molecule_chembl_id", "inchikey14"]]]).drop_duplicates("molecule_chembl_id")
        df = df.merge(ik, on="molecule_chembl_id", how="left")
        conflict = df.groupby("inchikey14")["y"].nunique()
        n_conflict = int((conflict > 1).sum())

        # 统一骨架五折
        df = df.reset_index(drop=True)
        df["fold"] = -1
        sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
        for k, (_, te) in enumerate(sgkf.split(df, df["y"], groups=df["scaffold"])):
            df.loc[te, "fold"] = k

        # 跨折检查
        cross_scaf = int((df.groupby("scaffold")["fold"].nunique() > 1).sum())
        cross_ik = int((df.groupby("inchikey14")["fold"].nunique() > 1).sum())

        d = df.copy()
        d["target"] = target
        split_rows.append(d[["target", "molecule_chembl_id", "y", "scaffold", "inchikey14", "fold"]])

        X = np.vstack([ecfp[m] for m in df["molecule_chembl_id"]]).astype(np.float32)
        y = df["y"].to_numpy()
        folds = df["fold"].to_numpy()
        oof = np.full(len(df), np.nan)
        for k in range(5):
            tr, te = folds != k, folds == k
            if te.sum() == 0 or len(np.unique(y[tr])) < 2:
                continue
            m = RandomForestClassifier(**params)
            m.fit(X[tr], y[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]
        mask = ~np.isnan(oof)
        mt = metrics(y[mask], oof[mask])
        prev = float(y[mask].mean())

        all_rows.append({"target": target, "benchmark": "hard_unified_folds",
                         "n_pos": int(y[mask].sum()), "n_neg": int((1 - y[mask]).sum()),
                         "prevalence": prev, "pr_auc": mt["pr_auc"], "pr_sd": np.nan,
                         "lift": mt["pr_auc"] / prev, "roc_auc": mt["roc_auc"],
                         "ef1": mt["ef1"], "p@10": mt["p@10"]})

        easy_row = perf[(perf.target == target) & (perf.model == "rf")
                        & (perf.features == "ecfp4")].iloc[0]
        all_rows.append({"target": target, "benchmark": "easy",
                         "n_pos": int(easy_row.get("n") or 0) or int((splits[(splits.target == target) & (splits.y == 1)]).shape[0]),
                         "n_neg": int((splits[(splits.target == target) & (splits.y == 0)]).shape[0]),
                         "prevalence": np.nan, "pr_auc": easy_row["pr_auc"], "pr_sd": np.nan,
                         "lift": np.nan, "roc_auc": easy_row["roc_auc"],
                         "ef1": easy_row["ef1"], "p@10": easy_row["p@10"]})

        print(f"{target}: hard_unified PR={mt['pr_auc']:.4f} ROC={mt['roc_auc']:.4f} "
              f"prev={prev:.3f} | 跨折 scaffold={cross_scaf} ik={cross_ik} 冲突={n_conflict} "
              f"剔除跨靶点={dropped}")
        log += [f"### {target}",
                f"- 合并集合：{len(df)}（正 {int(y.sum())} / 负 {int((1 - y).sum())}），"
                f"骨架数 {df['scaffold'].nunique()}",
                f"- 统一 StratifiedGroupKFold(5, seed={SEED}) 后跨折检查："
                f"scaffold 跨折 {cross_scaf}，InChIKey14 跨折 {cross_ik}，标签冲突结构 {n_conflict}",
                f"- 剔除跨靶点分子（骨架/结构冲突）{dropped}",
                f"- hard_unified：PR-AUC={mt['pr_auc']:.4f}，ROC-AUC={mt['roc_auc']:.4f}，"
                f"lift={mt['pr_auc'] / prev:.2f}",
                ""]

    splits_out = pd.concat(split_rows, ignore_index=True)
    splits_out.to_csv(PROC / "hard_benchmark_splits.csv", index=False, encoding="utf-8-sig")

    out_path = TAB / "ligand_model_hard_benchmark.csv"
    legacy = TAB / "ligand_model_hard_benchmark_legacy.csv"
    if not legacy.exists():
        shutil.copy(out_path, legacy)
    pd.DataFrame(all_rows).to_csv(out_path, index=False, encoding="utf-8-sig")

    log += ["- 原切分（各靶点独立折号）保留在 data/processed/scaffold_splits.csv；"
            "新版切分写入 data/processed/hard_benchmark_splits.csv",
            "- 其他靶点配体仅为“假定困难负例”，不等于该靶点实测非活性分子",
            "- prevalence 变化会影响 AP，prevalence-matched 下采样仅作辅助描述", ""]
    with open(LOGS / "p7c_hard_benchmark.md", "a", encoding="utf-8") as f:
        f.write("\n".join(log))

    print(f"[18b] 完成，耗时 {(time.time() - T0) / 60:.1f} 分钟")


if __name__ == "__main__":
    main()