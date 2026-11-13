# -*- coding: utf-8 -*-
"""P7b 诊断：泄漏体检、标签置换、描述符偏置、跨靶点硬负样本、候选新颖性、随机基线重标定。

输入：results/tables/pool_actives.csv、pool_decoys.csv、top_candidates.csv、
      data/interim/actives_raw.csv、decoys_raw.csv、data/processed/scaffold_splits.csv、
      ligand_features_ecfp4.npz、ligand_features_desc.csv、ligand_oof_predictions.csv、
      ligand_model_perf.csv、ablation_summary.csv
输出：results/tables/*_diagnostics.csv、descriptor_auc.csv、cross_target_negatives.csv、
      top_candidates.csv（新增新颖性两列）、logs/p7b_diagnostics.md
运行：C:/Anaconda/Anaconda3/python.exe src/17_diagnostics.py

约束：不改 scaffold_splits.csv、不改主模型配置、不覆盖 ligand_model_perf.csv。
若标签置换 PR-AUC > 0.30 立即停止后续步骤。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import (DESC_COLS, LOGS, MAIN_TARGETS, SEED, TAB,  # noqa: E402
                          load_features, load_splits)

RDLogger.DisableLog("rdApp.*")
SEED = 42


def ik14(smiles):
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) and smiles.strip() else None
    if m is None:
        return None
    try:
        return Chem.MolToInchiKey(m)[:14]
    except Exception:
        return None


def metrics(y, p):
    return float(roc_auc_score(y, p)), float(average_precision_score(y, p))


# ---------------- 一、集合与切分体检 ----------------
def leak_check():
    rows = []
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    raw_a = pd.read_csv(Path("data/interim/actives_raw.csv"), encoding="utf-8-sig")
    raw_d = pd.read_csv(Path("data/interim/decoys_raw.csv"), encoding="utf-8-sig")

    for t in act["target"].unique():
        a_post = set(act.loc[act.target == t, "inchikey14"])
        d_post = set(dec.loc[dec.target == t, "inchikey14"])
        a_pre = {ik14(s) for s in raw_a.loc[raw_a.target == t, "smiles"]} - {None}
        d_pre = {ik14(s) for s in raw_d.loc[raw_d.target == t, "smiles"]} - {None}
        rows.append({"check": "actives_vs_decoys_overlap", "target": t,
                     "pre_standardization": len(a_pre & d_pre),
                     "post_standardization": len(a_post & d_post), "detail": ""})
    return pd.DataFrame(rows)


def split_check():
    sp = pd.read_csv(Path("data/processed/scaffold_splits.csv"), encoding="utf-8-sig")
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    ik = pd.concat([act[["molecule_chembl_id", "inchikey14"]],
                    dec[["molecule_chembl_id", "inchikey14"]]]).drop_duplicates("molecule_chembl_id")
    sp = sp.merge(ik, on="molecule_chembl_id", how="left")

    rows = []
    for t in sp["target"].unique():
        sub = sp[sp.target == t]
        for key, name in [("group_key", "scaffold/group_key"), ("inchikey14", "inchikey14")]:
            g = sub.groupby(key)["fold"].nunique()
            n_cross = int((g > 1).sum())
            rows.append({"check": f"cross_fold_{name}", "target": t,
                         "n_cross_fold_keys": n_cross,
                         "detail": "" if n_cross == 0 else "分组键失效，需排查"})
    return pd.DataFrame(rows)


# ---------------- 二、标签置换 ----------------
def label_permutation(splits, ecfp):
    rows = []
    for t in MAIN_TARGETS:
        df = splits[splits.target == t].reset_index(drop=True)
        y = df["y"].to_numpy()
        folds = df["fold"].to_numpy()
        X = np.vstack([ecfp[m] for m in df["molecule_chembl_id"]]).astype(np.float32)
        rng = np.random.default_rng(SEED)
        oof = np.full(len(df), np.nan)
        for k in sorted(np.unique(folds)):
            tr, te = folds != k, folds == k
            y_perm = y[tr].copy()
            rng.shuffle(y_perm)
            if len(np.unique(y_perm)) < 2:
                continue
            m = RandomForestClassifier(n_estimators=500, class_weight="balanced_subsample",
                                       random_state=SEED, n_jobs=-1)
            m.fit(X[tr], y_perm)
            oof[te] = m.predict_proba(X[te])[:, 1]
        mask = ~np.isnan(oof)
        roc, pr = metrics(y[mask], oof[mask])
        rows.append({"target": t, "roc_auc": roc, "pr_auc": pr,
                     "prevalence": float(y.mean()), "n": int(len(df))})
        print(f"[permutation] {t}: ROC={roc:.4f} PR={pr:.4f} (prevalence={y.mean():.3f})")
    return pd.DataFrame(rows)


# ---------------- 三、描述符 AUC ----------------
def descriptor_auc(splits, desc):
    rows = []
    for t in MAIN_TARGETS:
        df = splits[splits.target == t]
        y = df["y"].to_numpy()
        D = desc.loc[df["molecule_chembl_id"].tolist(), DESC_COLS].to_numpy()
        for j, c in enumerate(DESC_COLS):
            auc = roc_auc_score(y, D[:, j])
            rows.append({"target": t, "descriptor": c, "auc": auc,
                         "auc_flipped": max(auc, 1 - auc)})
    return pd.DataFrame(rows)


# ---------------- 四、跨靶点硬负样本 ----------------
def cross_target(splits, ecfp, perf, oof_pred):
    rows = []
    # 每个靶点的 Youden 阈值来自 12 号 OOF
    cut = {}
    for t in MAIN_TARGETS:
        sub = oof_pred[(oof_pred.target == t) & (oof_pred.model == "rf")]
        from sklearn.metrics import roc_curve
        fpr, tpr, thr = roc_curve(sub["y"], sub["prob"])
        cut[t] = float(thr[int(np.argmax(tpr - fpr))])

    def fit_target(t):
        df = splits[splits.target == t].reset_index(drop=True)
        X = np.vstack([ecfp[m] for m in df["molecule_chembl_id"]]).astype(np.float32)
        y = df["y"].to_numpy()
        m = RandomForestClassifier(n_estimators=500, class_weight="balanced_subsample",
                                   random_state=SEED, n_jobs=-1)
        m.fit(X, y)
        return m, df

    for model_t, donor_ts in [("PDE4B", ["PFKFB3", "EPAS1"]), ("PFKFB3", ["PDE4B"])]:
        model, df_m = fit_target(model_t)
        own_ids = set(df_m.loc[df_m.y == 1, "molecule_chembl_id"])
        Xown = np.vstack([ecfp[m] for m in df_m["molecule_chembl_id"]]).astype(np.float32)
        p_own = model.predict_proba(Xown)[:, 1]

        for donor in donor_ts:
            dd = splits[(splits.target == donor) & (splits.y == 1)]
            # 扣除与模型靶点训练集 InChIKey14 重叠者
            act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
            own_ik = set(act[act.target == model_t]["inchikey14"])
            donor_ik = act[act.target == donor].set_index("molecule_chembl_id")["inchikey14"]
            keep = [m for m in dd["molecule_chembl_id"] if donor_ik.get(m) not in own_ik]
            if not keep:
                rows.append({"model_target": model_t, "donor_target": donor,
                             "n_negatives": 0, "frac_pred_pos": np.nan,
                             "median_percentile": np.nan})
                continue
            Xn = np.vstack([ecfp[m] for m in keep]).astype(np.float32)
            p = model.predict_proba(Xn)[:, 1]
            pct = [(p_own <= v).mean() * 100 for v in p]
            rows.append({"model_target": model_t, "donor_target": donor,
                         "n_negatives": len(keep),
                         "frac_pred_pos": float((p >= cut[model_t]).mean()),
                         "median_percentile": float(np.median(pct)),
                         "cutoff": round(cut[model_t], 4)})
            print(f"[cross] {model_t} 模型 <- {donor} 正样本 n={len(keep)}: "
                  f"frac_pred_pos={rows[-1]['frac_pred_pos']:.3f}, "
                  f"median_pct={rows[-1]['median_percentile']:.1f}")
    return pd.DataFrame(rows)


# ---------------- 五、候选新颖性 ----------------
def novelty(ecfp, splits):
    tc = pd.read_csv(TAB / "top_candidates.csv", encoding="utf-8-sig")
    out, counts = [], {}
    for t in MAIN_TARGETS:
        train_ids = splits[(splits.target == t) & (splits.y == 1)]["molecule_chembl_id"].tolist()
        train_fps = [ecfp[m] for m in train_ids]
        n_novel = 0
        for r in tc[tc.target == t].itertuples(index=False):
            v = ecfp[r.molecule_chembl_id].astype(np.float32)
            A = np.vstack(train_fps).astype(np.float32)
            inter = A @ v
            union = A.sum(axis=1) + v.sum() - inter
            best = float((inter / np.maximum(union, 1)).max())
            flag = "analog" if best > 0.8 else ("similar" if best >= 0.4 else "novel")
            if flag == "novel":
                n_novel += 1
            out.append({"target": t, "rank": r.rank, "molecule_chembl_id": r.molecule_chembl_id,
                        "pref_name": r.pref_name, "prob": r.prob,
                        "max_tanimoto_to_train_actives": round(best, 4), "novelty_flag": flag})
        counts[t] = n_novel
    res = pd.DataFrame(out)
    res.to_csv(TAB / "top_candidates_novelty_diagnostics.csv", index=False, encoding="utf-8-sig")

    tc2 = tc.merge(res[["target", "rank", "max_tanimoto_to_train_actives", "novelty_flag"]],
                   on=["target", "rank"], how="left")
    tc2.to_csv(TAB / "top_candidates.csv", index=False, encoding="utf-8-sig")
    return res, counts


# ---------------- 六、随机基线重标定 ----------------
def random_baseline(splits):
    from ligand_utils import ef_at, precision_at
    rows = []
    for t in MAIN_TARGETS:
        df = splits[splits.target == t]
        y = df["y"].to_numpy()
        n = len(y)
        rng = np.random.default_rng(SEED)
        prs, ef1s = [], []
        for _ in range(100):
            p = rng.random(n)
            prs.append(average_precision_score(y, p))
            ef1s.append(ef_at(y, p, 0.01))
        rows.append({"target": t, "n": n, "prevalence": float(y.mean()),
                     "pr_mean": float(np.mean(prs)), "pr_sd": float(np.std(prs)),
                     "ef1_mean": float(np.mean(ef1s)), "ef1_sd": float(np.std(ef1s))})
    return pd.DataFrame(rows)


def main():
    ecfp, desc = load_features()
    splits = load_splits()
    perf = pd.read_csv(TAB / "ligand_model_perf.csv", encoding="utf-8-sig")
    oof_pred = pd.read_csv(TAB / "ligand_oof_predictions.csv", encoding="utf-8-sig")
    ab = pd.read_csv(TAB / "ablation_summary.csv", encoding="utf-8-sig")

    logs = ["# P7b 诊断日志", ""]

    # 一
    lc = leak_check()
    sc = split_check()
    checks = pd.concat([lc, sc], ignore_index=True)
    checks.to_csv(TAB / "leak_check_diagnostics.csv", index=False, encoding="utf-8-sig")
    print(checks.to_string(index=False))
    logs += ["## 一、集合与切分体检", "", checks.to_string(index=False), ""]

    # 二
    perm = label_permutation(splits, ecfp)
    perm.to_csv(TAB / "label_permutation_diagnostics.csv", index=False, encoding="utf-8-sig")
    logs += ["## 二、标签置换对照", "", perm.round(4).to_string(index=False), ""]
    if (perm["pr_auc"] > 0.30).any():
        logs.append("**标签置换 PR-AUC > 0.30，判定存在特征/样本对齐问题，后续步骤已停止。**")
        (LOGS / "p7b_diagnostics.md").write_text("\n".join(logs), encoding="utf-8")
        print("[STOP] 置换 PR-AUC > 0.30，停止后续诊断")
        return

    # 三
    logs += ["## 三、消融全表", "", ab.round(4).to_string(index=False), ""]
    da = descriptor_auc(splits, desc)
    da.to_csv(TAB / "descriptor_auc.csv", index=False, encoding="utf-8-sig")
    best = da.sort_values("auc_flipped", ascending=False).head(5)
    logs += ["### 单维描述符 AUC（auc_flipped = 方向归一后的可分性）", "",
             best.round(4).to_string(index=False), ""]

    # 四
    ct = cross_target(splits, ecfp, perf, oof_pred)
    ct.to_csv(TAB / "cross_target_negatives.csv", index=False, encoding="utf-8-sig")
    logs += ["## 四、跨靶点硬负样本", "", ct.round(4).to_string(index=False), ""]

    # 五
    nov, counts = novelty(ecfp, splits)
    logs += ["## 五、Top 候选新颖性", "",
             f"- novel 个数：{counts}", ""]

    # 六
    rb = random_baseline(splits)
    rb.to_csv(TAB / "baseline_diagnostics.csv", index=False, encoding="utf-8-sig")
    logs += ["## 六、随机基线（100 次置换 mean±SD）", "", rb.round(4).to_string(index=False), ""]

    (LOGS / "p7b_diagnostics.md").write_text("\n".join(logs), encoding="utf-8")
    print("[17_diagnostics] 完成")


if __name__ == "__main__":
    main()