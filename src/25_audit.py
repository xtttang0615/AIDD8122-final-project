# -*- coding: utf-8 -*-
"""25 审查脚本：修复前的证据核查（A 分子身份对应 / C 与全部训练数据重叠）。

输出：results/tables/audit_mi_identity.csv、audit_overlap.csv、logs/revision_audit.md
运行：C:/Anaconda/Anaconda3/python.exe src/25_audit.py
只读核查，不修改任何既有结果表。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, MAIN_TARGETS, ROOT, TAB  # noqa: E402
from pool_utils import descriptors, prepare_smiles  # noqa: E402


def audit_mi_identity():
    """复现 14 号脚本的特征过滤顺序，判断 gold.iloc[:len(G_fp)] 是否发生身份错配。"""
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")
    gold = gold[gold["smiles"].notna() & (gold["smiles"].astype(str).str.strip() != "")]
    ok_flags, reasons = [], []
    for s in gold["smiles"]:
        mol, reason = prepare_smiles(s)
        ok_flags.append(mol is not None)
        reasons.append(reason or "")
    n_ok = int(sum(ok_flags))
    first_n_ids = set(gold["molecule_chembl_id"].iloc[:n_ok])
    ok_ids = set(gold.loc[pd.Series(ok_flags, index=gold.index), "molecule_chembl_id"])

    mismatch = first_n_ids != ok_ids
    failed = gold.loc[~pd.Series(ok_flags, index=gold.index), ["molecule_chembl_id", "pref_name"]]
    failed = failed.assign(reason=[r for r, ok in zip(reasons, ok_flags) if not ok])

    rows = pd.DataFrame([{
        "n_with_smiles": len(gold),
        "n_feature_ok": n_ok,
        "n_failed": len(failed),
        "identity_mismatch_by_truncation": mismatch,
        "n_ids_only_in_truncation": len(first_n_ids - ok_ids),
        "n_ids_only_in_success": len(ok_ids - first_n_ids),
        "failed_ids": ";".join(failed["molecule_chembl_id"].tolist()),
        "failed_names": ";".join(failed["pref_name"].astype(str).tolist()),
    }])
    return rows, failed


def audit_overlap():
    """候选池/评价集 与 每个靶点训练正例、训练负例（decoys）的结构键重叠。"""
    pool = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")
    scores = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")

    rows = []
    for name, df in [("approved_pool", pool), ("approved_scores", scores)]:
        for t in MAIN_TARGETS + ["PDK4"]:
            pos = set(act.loc[act.target == t, "inchikey14"])
            neg = set(dec.loc[dec.target == t, "inchikey14"])
            ik = set(df["inchikey14"])
            rows.append({
                "set": name, "target": t, "n_set": len(ik),
                "n_overlap_train_pos": len(ik & pos),
                "n_overlap_train_neg": len(ik & neg),
                "n_overlap_all_train": len(ik & (pos | neg)),
                "overlap_ids": ";".join(sorted(ik & (pos | neg))[:20]),
            })
    # MI 评价集
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")
    g_ik = set()
    for s in gold["smiles"].dropna():
        mol, _ = prepare_smiles(s)
        if mol is not None:
            g_ik.add(descriptors(mol)["inchikey14"])
    for t in MAIN_TARGETS + ["PDK4"]:
        pos = set(act.loc[act.target == t, "inchikey14"])
        neg = set(dec.loc[dec.target == t, "inchikey14"])
        rows.append({"set": "mi_gold", "target": t, "n_set": len(g_ik),
                     "n_overlap_train_pos": len(g_ik & pos),
                     "n_overlap_train_neg": len(g_ik & neg),
                     "n_overlap_all_train": len(g_ik & (pos | neg)),
                     "overlap_ids": ";".join(sorted(g_ik & (pos | neg))[:20])})
    return pd.DataFrame(rows)


def main():
    ident, failed = audit_mi_identity()
    ident.to_csv(TAB / "audit_mi_identity.csv", index=False, encoding="utf-8-sig")
    failed.to_csv(TAB / "audit_mi_identity_failed.csv", index=False, encoding="utf-8-sig")

    ov = audit_overlap()
    ov.to_csv(TAB / "audit_overlap.csv", index=False, encoding="utf-8-sig")

    print("=== A 分子身份核查 ===")
    print(ident.T.to_string())
    print("\n=== C 重叠核查 ===")
    print(ov.to_string(index=False))

    log = ["# 修订前核查记录（P9-audit）", "",
           "## A. 14 号脚本分子身份对应", "",
           ident.T.to_string(), "",
           "失败分子：", failed.to_string(index=False), "",
           "## C. 候选/评价集与训练数据重叠（InChIKey14）", "",
           ov.to_string(index=False), ""]
    (LOGS / "revision_audit.md").write_text("\n".join(log), encoding="utf-8")


if __name__ == "__main__":
    main()