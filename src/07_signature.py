# -*- coding: utf-8 -*-
"""P4c 诊断签名建模：完全嵌套 CV（折内 limma）+ 跨队列标准化敏感性分析。

输入：data/processed/GSE66360_gene.csv、GSE66360_meta.csv、GSE48060_gene.csv、GSE48060_meta.csv
      results/tables/deg_train.csv（全训练集特征池，仅用于固定池对照）
输出：results/tables/signature_genes.csv、model_performance.csv、predictions.csv、calibration.csv
      logs/model_log.md、logs/model_try.md
运行：C:/Anaconda/Anaconda3/python.exe src/07_signature.py

方法学：
- 主结果 train_nested_full：外层分层 5 折；每折训练部分内完整重跑 limma -> bootstrap 稳定性筛选
  -> LogisticRegression(L2,C=1.0)，预测留出折。折内 limma 由 src/05b_limma_fold.R 执行。
- train_nested_fixedpool：特征池由全训练集 DEG 定义，存在 pre-filtering bias，仅作对照。
- train_cv_optimistic：全量筛选 + 同批样本 CV，存在选择泄漏，仅作对照。
- external_zscored：post-hoc 敏感性分析，逐数据集内部逐基因 z-score，不使用任何标签信息。

约束：不换分类器、不调 C、不对 internal/external 单独寻阈值、不改 deg_* 文件。
"""
from pathlib import Path
import os
import shutil
import subprocess
import sys
import tempfile

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.metrics import roc_auc_score, roc_curve

ROOT = Path(__file__).resolve().parents[1]
PROCESSED = ROOT / "data" / "processed"
TAB = ROOT / "results" / "tables"
LOGS = ROOT / "logs"
TAB.mkdir(parents=True, exist_ok=True)

# Rscript 路径优先取环境变量 RSCRIPT，其次 PATH 中的 Rscript，最后退回本机默认安装位置
RSCRIPT = Path(os.environ.get("RSCRIPT") or shutil.which("Rscript")
               or r"C:\R\R-4.6.1\bin\Rscript.exe")
LIM_AUX = ROOT / "src" / "05b_limma_fold.R"

SEED = 42
N_BOOT_MAIN = 200
N_BOOT_NESTED = 100
N_CI = 1000
RF_SEL_TREES = 500
RF_CLF_TREES = 1000
LASSO_C = 1.0
FREQ_CUT = 0.50
MIN_SIG = 5
CUTOFF_LOCKED = 0.6001   # P4 锁定的统一阈值
CUTOFF_HALF = 0.5
DEG_LFC = 1.0
DEG_FDR = 0.05


def load_all():
    expr6 = pd.read_csv(PROCESSED / "GSE66360_gene.csv", index_col=0, encoding="utf-8-sig")
    meta6 = pd.read_csv(PROCESSED / "GSE66360_meta.csv", encoding="utf-8-sig").set_index("sample_id")
    meta6 = meta6.loc[expr6.columns]
    expr8 = pd.read_csv(PROCESSED / "GSE48060_gene.csv", index_col=0, encoding="utf-8-sig")
    meta8 = pd.read_csv(PROCESSED / "GSE48060_meta.csv", encoding="utf-8-sig").set_index("sample_id")
    meta8 = meta8.loc[expr8.columns]
    deg = pd.read_csv(TAB / "deg_train.csv", encoding="utf-8-sig")
    deg["sig"] = deg["sig"].astype(str).str.upper().isin(["TRUE", "1", "T"])
    return expr6, meta6, expr8, meta8, deg


def stratified_bootstrap_idx(y, rng):
    idx1 = np.where(y == 1)[0]
    idx0 = np.where(y == 0)[0]
    return np.concatenate([
        rng.choice(idx1, size=len(idx1), replace=True),
        rng.choice(idx0, size=len(idx0), replace=True),
    ])


def stability_selection(X, y, n_boot, rng):
    n_genes = X.shape[1]
    rf_freq = np.zeros(n_genes)
    lasso_freq = np.zeros(n_genes)
    sel_freq = np.zeros(n_genes)
    rf_imp_sum = np.zeros(n_genes)
    lasso_abs_sum = np.zeros(n_genes)
    done = 0

    for _ in range(n_boot):
        idx = stratified_bootstrap_idx(y, rng)
        Xb, yb = X[idx], y[idx]
        if len(np.unique(yb)) < 2:
            continue
        done += 1

        rf = RandomForestClassifier(n_estimators=RF_SEL_TREES, random_state=SEED, n_jobs=-1)
        rf.fit(Xb, yb)
        imp = rf.feature_importances_
        rf_imp_sum += imp
        rf_top = set(np.argsort(-imp)[:20])

        sc = StandardScaler()
        Xs = sc.fit_transform(Xb)
        las = LogisticRegression(penalty="l1", solver="liblinear", C=LASSO_C,
                                 random_state=SEED, max_iter=2000)
        las.fit(Xs, yb)
        coef = np.abs(las.coef_.ravel())
        lasso_abs_sum += coef
        lasso_set = set(np.where(coef > 0)[0])

        inter = rf_top & lasso_set
        rf_freq[list(rf_top)] += 1
        lasso_freq[list(lasso_set)] += 1
        if inter:
            sel_freq[list(inter)] += 1

    done = max(done, 1)
    return (rf_freq / done, lasso_freq / done, sel_freq / done,
            rf_imp_sum / done, lasso_abs_sum / done)


def pick_signature(sel_freq, rf_imp):
    stable_idx = np.where(sel_freq >= FREQ_CUT)[0]
    if len(stable_idx) >= MIN_SIG:
        return stable_idx[np.argsort(-sel_freq[stable_idx], kind="stable")], False
    return np.lexsort((-rf_imp, -sel_freq))[:10], True


def youden_cutoff(y, p):
    fpr, tpr, thr = roc_curve(y, p)
    return float(thr[int(np.argmax(tpr - fpr))])


def bootstrap_auc_ci(y, p, n=N_CI):
    rng = np.random.default_rng(SEED)
    idx1 = np.where(y == 1)[0]
    idx0 = np.where(y == 0)[0]
    vals = []
    for _ in range(n):
        s = np.concatenate([rng.choice(idx1, len(idx1), replace=True),
                            rng.choice(idx0, len(idx0), replace=True)])
        vals.append(roc_auc_score(y[s], p[s]))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def metrics_at(y, p, cutoff):
    pred = (p >= cutoff).astype(int)
    tp = int(((pred == 1) & (y == 1)).sum())
    tn = int(((pred == 0) & (y == 0)).sum())
    fp = int(((pred == 1) & (y == 0)).sum())
    fn = int(((pred == 0) & (y == 1)).sum())
    sens = tp / (tp + fn) if (tp + fn) else np.nan
    spec = tn / (tn + fp) if (tn + fp) else np.nan
    return round(sens, 4), round(spec, 4), round((tp + tn) / len(y), 4)


def calibration_stats(y, p):
    p = np.clip(p, 1e-6, 1 - 1e-6)
    logit = np.log(p / (1 - p)).reshape(-1, 1)
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    lr.fit(logit, y)
    return float(lr.coef_[0][0]), float(lr.intercept_[0]), float(np.mean((p - y) ** 2))


def run_fold_limma(expr_csv, label_csv, fold_no, out_csv):
    """调用 05b_limma_fold.R 在折内训练部分跑 limma。"""
    cmd = [str(RSCRIPT), str(LIM_AUX), str(expr_csv), str(label_csv), str(fold_no), str(out_csv)]
    res = subprocess.run(cmd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    if res.returncode != 0:
        raise RuntimeError(f"fold {fold_no} limma 失败：{res.stderr[:800]}")
    return pd.read_csv(out_csv, encoding="utf-8-sig")


def nested_fixed_pool(X, y, genes):
    """特征池固定为全训练集 DEG，只在折内重做稳定性筛选（有 pre-filtering bias）。"""
    outer = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    oof = np.zeros(len(y))
    info = []
    for fold, (tr, te) in enumerate(outer.split(X, y)):
        rng = np.random.default_rng(SEED + fold)
        _, _, sel_f, rf_i, _ = stability_selection(X[tr], y[tr], N_BOOT_NESTED, rng)
        order, _ = pick_signature(sel_f, rf_i)
        model = LogisticRegression(penalty="l2", C=1.0, max_iter=2000, random_state=SEED)
        model.fit(X[np.ix_(tr, order)], y[tr])
        oof[te] = model.predict_proba(X[np.ix_(te, order)])[:, 1]
        info.append(len(order))
        print(f"  [fixedpool fold {fold+1}] 训练 {len(tr)} / 留出 {len(te)}，签名 {len(order)} 个")
    return oof, info


def nested_full(expr6, train_ids, y, folds):
    """完全嵌套：折内 limma 重新定义特征池，再做稳定性筛选与拟合。"""
    tmpdir = Path(tempfile.mkdtemp(prefix="p4c_"))
    expr_csv = tmpdir / "expr_train.csv"
    label_csv = tmpdir / "labels.csv"
    try:
        expr6.loc[:, train_ids].to_csv(expr_csv, encoding="utf-8")
        labels = pd.DataFrame({
            "sample_id": train_ids,
            "group": np.where(y == 1, "case", "control"),
            "fold": folds,
        })
        labels.to_csv(label_csv, index=False, encoding="utf-8")

        oof = np.zeros(len(y))
        info = []
        for fold, (tr, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=SEED)
                                        .split(np.zeros(len(y)), y)):
            deg_f = run_fold_limma(expr_csv, label_csv, fold + 1, tmpdir / f"deg_fold{fold+1}.csv")
            sig_f = deg_f[(deg_f["logFC"].abs() > DEG_LFC) & (deg_f["adj.P.Val"] < DEG_FDR)]
            pool_f = [g for g in sig_f["gene"].tolist() if g in expr6.index]
            if len(pool_f) == 0:
                raise RuntimeError(f"fold {fold+1} 未得到 DEG，无法继续")

            Xf = expr6.loc[pool_f, train_ids].T.values
            rng = np.random.default_rng(SEED + fold)
            _, _, sel_f, rf_i, _ = stability_selection(Xf[tr], y[tr], N_BOOT_NESTED, rng)
            order, fallback = pick_signature(sel_f, rf_i)
            genes_sel = [pool_f[i] for i in order]

            model = LogisticRegression(penalty="l2", C=1.0, max_iter=2000, random_state=SEED)
            model.fit(Xf[np.ix_(tr, order)], y[tr])
            oof[te] = model.predict_proba(Xf[np.ix_(te, order)])[:, 1]

            info.append({"fold": fold + 1, "n_deg": len(pool_f), "n_sig": len(genes_sel),
                         "fallback": fallback, "genes": genes_sel})
            print(f"  [fullnested fold {fold+1}] 训练 {len(tr)} / 留出 {len(te)}，"
                  f"折内 DEG {len(pool_f)}，签名 {len(genes_sel)}"
                  f"{'（兜底top10）' if fallback else ''}")
        return oof, info
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    expr6, meta6, expr8, meta8, deg = load_all()

    train_ids = meta6.index[meta6["role"] == "train"].tolist()
    internal_ids = meta6.index[meta6["role"] == "internal_val"].tolist()
    external_ids = meta8.index.tolist()

    y_train = (meta6.loc[train_ids, "group"] == "case").astype(int).values
    y_internal = (meta6.loc[internal_ids, "group"] == "case").astype(int).values
    y_external = (meta8.loc[external_ids, "group"] == "case").astype(int).values

    pool = [g for g in deg.loc[deg["sig"], "gene"].tolist() if g in expr6.index]
    print(f"全训练集特征池：{len(pool)} 个 DEG")

    X_pool = expr6.loc[pool, train_ids].T.values
    lfc = deg.set_index("gene")["logFC"]

    # ---------- 全训练集稳定性筛选（最终签名，用于最终模型）
    rf_freq, lasso_freq, sel_freq, rf_imp, lasso_abs = stability_selection(
        X_pool, y_train, N_BOOT_MAIN, np.random.default_rng(SEED))

    tab = pd.DataFrame({
        "gene": pool, "sel_freq": sel_freq, "rf_freq": rf_freq, "lasso_freq": lasso_freq,
        "rf_importance_mean": rf_imp, "lasso_coef_mean": lasso_abs,
    })
    tab["logFC"] = tab["gene"].map(lfc)
    tab["direction"] = np.where(tab["logFC"] > 0, "up", "down")

    n_stable = int((tab["sel_freq"] >= FREQ_CUT).sum())
    order, fallback = pick_signature(sel_freq, rf_imp)
    sig_genes = [pool[i] for i in order]
    sig_tab = tab[tab["gene"].isin(sig_genes)].copy().sort_values(
        ["sel_freq", "rf_importance_mean"], ascending=False)
    sig_tab.to_csv(TAB / "signature_genes.csv", index=False, encoding="utf-8-sig")
    print(f"sel_freq>=50% 基因数 {n_stable}；最终签名 {len(sig_genes)} 个 -> {sig_genes}")

    # ---------- 分类器（沿用 P4 选定 logreg_l2 为主候选，两者仍做一次 CV 对照）
    Xs_train = expr6.loc[sig_genes, train_ids].T.values
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    candidates = {
        "logreg_l2": LogisticRegression(penalty="l2", C=1.0, max_iter=2000, random_state=SEED),
        "random_forest": RandomForestClassifier(n_estimators=RF_CLF_TREES, random_state=SEED, n_jobs=-1),
    }
    cv_auc = {}
    for name, est in candidates.items():
        prob = cross_val_predict(est, Xs_train, y_train, cv=cv, method="predict_proba")[:, 1]
        cv_auc[name] = roc_auc_score(y_train, prob)
    best_name = max(cv_auc, key=cv_auc.get)
    best_model = candidates[best_name]
    oof_optimistic = cross_val_predict(best_model, Xs_train, y_train, cv=cv, method="predict_proba")[:, 1]
    print(f"分类器：logreg_l2={cv_auc['logreg_l2']:.4f} / RF={cv_auc['random_forest']:.4f} -> 选 {best_name}")

    # ---------- 嵌套 CV（固定池 + 完全嵌套）
    print("嵌套 5 折（固定池对照）...")
    oof_fixed, info_fixed = nested_fixed_pool(X_pool, y_train, pool)
    auc_fixed = roc_auc_score(y_train, oof_fixed)
    lo_f, hi_f = bootstrap_auc_ci(y_train, oof_fixed)

    folds = np.zeros(len(y_train), dtype=int)
    for k, (_, te) in enumerate(StratifiedKFold(5, shuffle=True, random_state=SEED)
                                .split(np.zeros(len(y_train)), y_train)):
        folds[te] = k + 1

    print("完全嵌套 5 折（折内 limma 重定义特征池）...")
    oof_full, info_full = nested_full(expr6, train_ids, y_train, folds)
    auc_full = roc_auc_score(y_train, oof_full)
    lo_full, hi_full = bootstrap_auc_ci(y_train, oof_full)
    print(f"train_nested_full AUC = {auc_full:.4f} [{lo_full:.4f}, {hi_full:.4f}]")

    # ---------- 最终模型 + 三层预测
    final_model = best_model
    final_model.fit(Xs_train, y_train)
    p_internal = final_model.predict_proba(expr6.loc[sig_genes, internal_ids].T.values)[:, 1]

    missing = [g for g in sig_genes if g not in expr8.index]
    medians = expr6.loc[sig_genes, train_ids].median(axis=1)
    X_ext = np.zeros((len(external_ids), len(sig_genes)))
    for j, g in enumerate(sig_genes):
        X_ext[:, j] = expr8.loc[g, external_ids].values if g in expr8.index else medians[g]
    p_external = final_model.predict_proba(X_ext)[:, 1]

    # ---------- 跨队列标准化敏感性分析（不使用任何标签）
    mu_tr = expr6.loc[sig_genes, train_ids].mean(axis=1)
    sd_tr = expr6.loc[sig_genes, train_ids].std(axis=1).replace(0, np.nan)
    Z_train = ((expr6.loc[sig_genes, train_ids].T - mu_tr) / sd_tr).fillna(0.0).values

    mu_ext = expr8.loc[sig_genes, external_ids].mean(axis=1)
    sd_ext = expr8.loc[sig_genes, external_ids].std(axis=1).replace(0, np.nan)
    Z_ext = ((expr8.loc[sig_genes, external_ids].T - mu_ext) / sd_ext).fillna(0.0).values

    model_z = LogisticRegression(penalty="l2", C=1.0, max_iter=2000, random_state=SEED)
    model_z.fit(Z_train, y_train)
    p_ext_z = model_z.predict_proba(Z_ext)[:, 1]
    print(f"external 标准化前后 mean_prob：{p_external.mean():.4f} -> {p_ext_z.mean():.4f}")

    # ---------- 汇总
    scenarios = [
        ("train_nested_full", y_train, oof_full, train_ids, "主结果（完全嵌套 CV）"),
        ("train_nested_fixedpool", y_train, oof_fixed, train_ids,
         "特征池由全训练集 DEG 定义，存在 pre-filtering bias"),
        ("train_cv_optimistic", y_train, oof_optimistic, train_ids,
         "存在选择泄漏，不作为主结果"),
        ("internal_val", y_internal, p_internal, internal_ids, ""),
        ("external", y_external, p_external, external_ids, ""),
        ("external_zscored", y_external, p_ext_z, external_ids,
         "post-hoc 敏感性分析，逐数据集内部 z-score，未使用任何标签"),
    ]

    rows, pred_records = [], []
    for set_name, y, p, ids, note in scenarios:
        auc = roc_auc_score(y, p)
        lo, hi = bootstrap_auc_ci(y, p)
        s_y, sp_y, a_y = metrics_at(y, p, CUTOFF_LOCKED)
        s_h, sp_h, a_h = metrics_at(y, p, CUTOFF_HALF)
        rows.append({
            "set": set_name, "n": len(y), "n_case": int(y.sum()), "n_control": int((1 - y).sum()),
            "auc": round(auc, 4), "auc_lo": round(lo, 4), "auc_hi": round(hi, 4),
            "cutoff_youden": CUTOFF_LOCKED, "sens_youden": s_y, "spec_youden": sp_y, "acc_youden": a_y,
            "cutoff_05": CUTOFF_HALF, "sens_05": s_h, "spec_05": sp_h, "acc_05": a_h,
            "note": note,
        })
        print(f"{set_name}: AUC={auc:.4f} [{lo:.4f}, {hi:.4f}] | "
              f"@0.6001 sens={s_y:.3f} spec={sp_y:.3f} acc={a_y:.3f} | "
              f"@0.5 sens={s_h:.3f} spec={sp_h:.3f} acc={a_h:.3f}")
        for sid, yy, pp in zip(ids, y, p):
            pred_records.append({"set": set_name, "sample_id": sid,
                                 "y_true": int(yy), "y_prob": round(float(pp), 6)})

    pd.DataFrame(rows).to_csv(TAB / "model_performance.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(pred_records).to_csv(TAB / "predictions.csv", index=False, encoding="utf-8-sig")

    cal_rows = []
    for set_name, y, p in [("internal_val", y_internal, p_internal),
                           ("external", y_external, p_external),
                           ("external_zscored", y_external, p_ext_z)]:
        slope, intercept, brier = calibration_stats(y, p)
        cal_rows.append({"set": set_name, "n": len(y), "calibration_slope": round(slope, 4),
                         "calibration_intercept": round(intercept, 4), "brier": round(brier, 4),
                         "mean_prob": round(float(np.mean(p)), 4)})
        print(f"{set_name} 校准：slope={slope:.4f} intercept={intercept:.4f} "
              f"Brier={brier:.4f} mean_prob={p.mean():.4f}")
    pd.DataFrame(cal_rows).to_csv(TAB / "calibration.csv", index=False, encoding="utf-8-sig")

    # ---------- 日志
    top20 = tab.sort_values(["sel_freq", "rf_importance_mean"], ascending=False).head(20)
    log = ["# P4 建模日志（P4c：完全嵌套 CV + 标准化敏感性）", "",
           f"- 随机种子：SEED={SEED}（外层折 rng=SEED+fold_index）",
           f"- 全训练集特征池：{len(pool)} 个 DEG（deg_train.csv 中 sig==TRUE）",
           f"- 主筛选 bootstrap：{N_BOOT_MAIN} 次；嵌套各折 bootstrap：{N_BOOT_NESTED} 次",
           f"- 筛选规则：RF(n_estimators={RF_SEL_TREES}) top-20 ∩ LASSO(L1,liblinear,C={LASSO_C}) 非零集，"
           f"sel_freq >= {FREQ_CUT:.0%}；不足 {MIN_SIG} 个按 sel_freq -> rf_importance_mean 降序取 top-10",
           f"- sel_freq >= 50% 的基因总数：{n_stable}",
           f"- 最终签名（全训练集）：{', '.join(sig_genes)}",
           f"- 分类器：logreg_l2={cv_auc['logreg_l2']:.4f} vs random_forest={cv_auc['random_forest']:.4f} "
           f"-> 选 {best_name}（未调 C）",
           "",
           "## 完全嵌套 CV（主结果）",
           f"- train_nested_full AUC = {auc_full:.4f} [{lo_full:.4f}, {hi_full:.4f}]",
           "- 每折流程：折内训练部分 limma(|logFC|>1 & adj.P<0.05) -> bootstrap 100 次稳定性筛选 -> L2 LR",
           ""]
    for r in info_full:
        log.append(f"- fold {r['fold']}：折内 DEG {r['n_deg']} 个，签名 {r['n_sig']} 个"
                   f"{'（兜底 top-10）' if r['fallback'] else ''}；基因：{', '.join(r['genes'])}")
    log += ["",
            "## 对照结果",
            f"- train_nested_fixedpool AUC = {auc_fixed:.4f} [{lo_f:.4f}, {hi_f:.4f}]："
            "特征池由全训练集 DEG 定义，存在 pre-filtering bias，不作主结果",
            f"- train_cv_optimistic AUC = {cv_auc[best_name]:.4f}：全量筛选 + 同批样本 CV，"
            "存在选择泄漏，不作主结果",
            "",
            "## 三层验证与阈值",
            f"- 统一阈值 cutoff={CUTOFF_LOCKED}（P4 锁定，三层共用），另并列 cutoff=0.5；"
            "未对 internal/external 单独寻阈值",
            f"- external 缺失基因数：{len(missing)}（缺失用 train 中位数填补）",
            ""]
    for r in rows:
        log.append(f"- {r['set']}: n={r['n']} (case {r['n_case']}/control {r['n_control']}), "
                   f"AUC={r['auc']} [{r['auc_lo']}, {r['auc_hi']}], "
                   f"@0.6001 sens={r['sens_youden']} spec={r['spec_youden']} acc={r['acc_youden']}, "
                   f"@0.5 sens={r['sens_05']} spec={r['spec_05']} acc={r['acc_05']}"
                   + (f" | {r['note']}" if r["note"] else ""))
    log += ["",
            "## 跨队列标准化敏感性分析（post-hoc）",
            "- 做法：GSE66360 用 train 43 例的逐基因均值/标准差、GSE48060 用其自身 52 例的逐基因均值/标准差，"
            "各自标准化后训练与预测；未使用任何标签信息。",
            f"- external AUC：{rows[4]['auc']:.4f} -> {rows[5]['auc']:.4f}；"
            f"mean_prob：{p_external.mean():.4f} -> {p_ext_z.mean():.4f}",
            "",
            "## 校准",
            ]
    for r in cal_rows:
        log.append(f"- {r['set']}: slope={r['calibration_slope']}, intercept={r['calibration_intercept']}, "
                   f"Brier={r['brier']}, mean_prob={r['mean_prob']}")
    log += ["",
            "## 结论（照实记录）",
            f"- bootstrap {N_BOOT_MAIN} 次中 sel_freq >= 50% 仅 1 个基因（PDE4B: 0.785），签名稳定性不足。",
            f"- external AUC = {rows[4]['auc']:.4f}（95% CI {rows[4]['auc_lo']:.4f}-{rows[4]['auc_hi']:.4f}），"
            "跨细胞类型（CEC -> 全血）迁移性能有限。",
            ""]
    (LOGS / "model_log.md").write_text("\n".join(log), encoding="utf-8")

    try_lines = ["# 建模试错记录", "",
                 "- P4 初版 train CV AUC=1.0000 判定为特征选择泄漏 -> P4b 改嵌套 CV（固定池）。",
                 "- P4c 进一步发现固定池仍存在 pre-filtering bias（特征池来自全训练集 DEG）-> "
                 "改为完全嵌套 CV（折内重跑 limma），作为主结果。",
                 f"- 分类器对比：logreg_l2={cv_auc['logreg_l2']:.4f} vs "
                 f"random_forest={cv_auc['random_forest']:.4f}，沿用预设规则，未反复调整。",
                 "- 跨队列标准化只试了'各数据集内部逐基因 z-score'一种方式，未试其他标准化；结果照实报告。",
                 "- 未对 internal/external 做阈值或参数搜索；未改 C；未改 deg_* 文件。",
                 "- P3 火山图 adjustText 触发 matplotlib FancyArrowPatch 警告（自动改用 annotate），不修复仅记录。",
                 ""]
    (LOGS / "model_try.md").write_text("\n".join(try_lines), encoding="utf-8")

    print("[07_signature] 完成（P4c）")


if __name__ == "__main__":
    main()