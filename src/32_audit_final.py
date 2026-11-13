# -*- coding: utf-8 -*-
"""32 提交前审核：核对报告关键数字与实际结果文件是否一致。只读，不改动任何结果。"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"


def n(p):
    return len(pd.read_csv(TAB / p, encoding="utf-8-sig"))


def main():
    meta = pd.read_csv(ROOT / "data/processed/GSE66360_meta.csv", encoding="utf-8-sig")
    print("role counts:", meta["role"].value_counts().to_dict())
    print("GSE48060 n:", len(pd.read_csv(ROOT / "data/processed/GSE48060_meta.csv", encoding="utf-8-sig")))

    g = pd.read_csv(ROOT / "data/processed/GSE66360_gene.csv", encoding="utf-8-sig", index_col=0)
    print("genes:", g.shape[0], "samples:", g.shape[1])

    deg = pd.read_csv(TAB / "deg_train.csv", encoding="utf-8-sig")
    print("DEG sig:", int(deg["sig"].sum()), "/ rows", len(deg))
    sig = pd.read_csv(TAB / "signature_genes.csv", encoding="utf-8-sig")
    print("signature genes:", len(sig), list(sig["gene"]))

    mp = pd.read_csv(TAB / "model_performance.csv", encoding="utf-8-sig")
    print(mp[["set", "auc", "auc_lo", "auc_hi"]].to_string(index=False))

    pa = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    pdec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    print("actives:", len(pa), pa["target"].value_counts().to_dict())
    print("decoys:", len(pdec), pdec["target"].value_counts().to_dict())

    ap = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")
    print("approved_pool rows:", len(ap))

    sc = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    print("scores rows:", len(sc), "mi_gold:", int(sc["is_mi_gold"].sum()))

    mi = pd.read_csv(TAB / "mi_enrichment_final.csv", encoding="utf-8-sig")
    print(mi.to_string(index=False))

    ctl = pd.read_csv(TAB / "supp_S3_control_molecules.csv", encoding="utf-8-sig")
    print(ctl.to_string(index=False))

    hb = pd.read_csv(TAB / "ligand_model_hard_benchmark.csv", encoding="utf-8-sig")
    print(hb[["target", "benchmark", "n_pos", "n_neg", "pr_auc", "roc_auc"]].to_string(index=False))

    ab = pd.read_csv(TAB / "ablation_summary.csv", encoding="utf-8-sig")
    print("ablation rows:", len(ab), "| no-op marked:", int((ab.get("note", pd.Series(dtype=str)).astype(str).str.contains("no-op|no_op", case=False)).sum()))


if __name__ == "__main__":
    main()