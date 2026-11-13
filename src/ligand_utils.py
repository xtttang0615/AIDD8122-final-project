# -*- coding: utf-8 -*-
"""P7 共用：特征加载、切分加载、模型工厂、OOF 指标。"""
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
PROC = ROOT / "data" / "processed"
LOGS = ROOT / "logs"

SEED = 42
NPZ = PROC / "ligand_features_ecfp4.npz"
DESC_CSV = PROC / "ligand_features_desc.csv"
SPLIT_CSV = PROC / "scaffold_splits.csv"

DESC_COLS = ["mw", "clogp", "hbd", "hba", "tpsa", "rotb", "rings", "arom_rings", "fsp3"]
MAIN_TARGETS = ["PDE4B", "PFKFB3", "EPAS1"]
ALL_TARGETS = MAIN_TARGETS + ["PDK4"]


def load_features():
    """返回 (ecfp4: dict id->np.uint8[2048], desc: DataFrame indexed by id)."""
    z = np.load(NPZ, allow_pickle=False)
    ids = z["ids"].astype(str)
    X = z["ecfp4"]
    ecfp = {m: X[i] for i, m in enumerate(ids)}
    desc = pd.read_csv(DESC_CSV, encoding="utf-8-sig").set_index("molecule_chembl_id")
    desc = desc.loc[[i for i in ids if i in desc.index]]
    return ecfp, desc


def load_splits():
    return pd.read_csv(SPLIT_CSV, encoding="utf-8-sig")


def load_pool(target):
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    a = act[act["target"] == target][["molecule_chembl_id"]].copy()
    a["y"] = 1
    d = dec[dec["target"] == target][["molecule_chembl_id"]].copy()
    d["y"] = 0
    return pd.concat([a, d], ignore_index=True)


def make_model(name, y_train=None):
    if name == "lr":
        return LogisticRegression(C=1.0, class_weight="balanced", max_iter=3000,
                                  random_state=SEED)
    if name == "rf":
        return RandomForestClassifier(n_estimators=500, class_weight="balanced_subsample",
                                      random_state=SEED, n_jobs=-1)
    if name == "xgb":
        from xgboost import XGBClassifier
        pos = max(int(np.sum(y_train == 1)), 1)
        neg = max(int(np.sum(y_train == 0)), 1)
        return XGBClassifier(n_estimators=300, max_depth=6, learning_rate=0.1,
                             scale_pos_weight=neg / pos, random_state=SEED,
                             eval_metric="logloss", n_jobs=-1, tree_method="hist")
    raise ValueError(name)


def ef_at(y, p, frac):
    """富集倍数：top frac 中的命中率 / 全体正例率。"""
    n = len(y)
    k = max(1, int(round(frac * n)))
    idx = np.argsort(-p)[:k]
    base = y.mean()
    if base <= 0:
        return np.nan
    return float((y[idx].sum() / k) / base)


def precision_at(y, p, k=10):
    idx = np.argsort(-p)[:k]
    return float(y[idx].sum() / min(k, len(y)))


def metrics(y, p):
    from sklearn.metrics import roc_auc_score, average_precision_score
    return {
        "roc_auc": float(roc_auc_score(y, p)),
        "pr_auc": float(average_precision_score(y, p)),
        "ef1": ef_at(y, p, 0.01),
        "ef5": ef_at(y, p, 0.05),
        "p@10": precision_at(y, p, 10),
    }


def bootstrap_ci(y, p, n_boot=1000, seed=SEED):
    from sklearn.metrics import roc_auc_score, average_precision_score
    rng = np.random.default_rng(seed)
    idx = np.arange(len(y))
    rocs, prs = [], []
    for _ in range(n_boot):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[s])) < 2:
            continue
        rocs.append(roc_auc_score(y[s], p[s]))
        prs.append(average_precision_score(y[s], p[s]))
    return (float(np.percentile(rocs, 2.5)), float(np.percentile(rocs, 97.5)),
            float(np.percentile(prs, 2.5)), float(np.percentile(prs, 97.5)))
