# -*- coding: utf-8 -*-
"""10b 拉取心肌梗死适应症药物作为评估金标准（仅用于评估，绝不进训练集）。

输入：ChEMBL drug_indication 接口
输出：results/tables/mi_gold_drugs.csv、logs/mi_drugs_log.md
运行：C:/Anaconda/Anaconda3/python.exe src/10b_fetch_mi_drugs.py

查询策略：先试 efo_id=EFO_0000612；若返回 0 条，改用 efo_term__contains="myocardial infarction"。
实际生效的查询方式写入日志与输出表。
"""
from pathlib import Path
import sys
from collections import Counter

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pool_utils import LOGS, ROOT, chembl_get, log_request_note  # noqa: E402

TAB = ROOT / "results" / "tables"
PAGE = 1000
CHUNK = 100


def fetch_indications(params, cache_prefix):
    rows, offset = [], 0
    while True:
        data, _ = chembl_get("/drug_indication.json",
                             {**params, "limit": PAGE, "offset": offset},
                             cache_name=f"{cache_prefix}_off{offset}.json")
        total = data["page_meta"]["total_count"]
        items = data.get("drug_indications") or data.get("drug_indication") or []
        rows.extend(items)
        offset += PAGE
        if not items or offset >= total:
            break
    return rows, total


def fetch_molecules(ids):
    out = {}
    ids = sorted(set(ids))
    for i in range(0, len(ids), CHUNK):
        chunk = ids[i:i + CHUNK]
        data, _ = chembl_get("/molecule.json",
                             {"molecule_chembl_id__in": ",".join(chunk), "limit": CHUNK,
                              "only": "molecule_chembl_id,pref_name,molecule_structures,max_phase"})
        for m in data.get("molecules") or []:
            st = m.get("molecule_structures") or {}
            out[m["molecule_chembl_id"]] = {
                "pref_name": m.get("pref_name"),
                "smiles": st.get("canonical_smiles"),
                "max_phase": m.get("max_phase"),
            }
    return out


def main():
    log_request_note("10b 开始拉取 MI 适应症药物")

    attempts = [
        ("efo_id=EFO_0000612", {"efo_id": "EFO_0000612"}, "di_efo"),
        ("efo_term__contains=myocardial infarction",
         {"efo_term__contains": "myocardial infarction"}, "di_term"),
    ]

    used, inds = None, []
    for label, params, prefix in attempts:
        rows, total = fetch_indications(params, prefix)
        print(f"[尝试] {label}: 命中 {total} 条")
        if rows:
            used, inds = label, rows
            break

    if not inds:
        raise RuntimeError("两种 drug_indication 查询方式均无结果，停止（不编造数据）")

    df = pd.DataFrame(inds)
    mols = fetch_molecules(df["molecule_chembl_id"].unique())

    df["pref_name"] = df["molecule_chembl_id"].map(lambda x: mols.get(x, {}).get("pref_name"))
    df["smiles"] = df["molecule_chembl_id"].map(lambda x: mols.get(x, {}).get("smiles"))
    df["max_phase"] = df["molecule_chembl_id"].map(lambda x: mols.get(x, {}).get("max_phase"))

    out = (df.groupby("molecule_chembl_id")
             .agg(pref_name=("pref_name", "first"), smiles=("smiles", "first"),
                  max_phase=("max_phase", "first"),
                  efo_terms=("efo_term", lambda s: ";".join(sorted(set(map(str, s))))),
                  n_indications=("efo_term", "size"))
             .reset_index())
    out.insert(0, "query_method", used)
    # ChEMBL 返回的 max_phase 是字符串（如 "4.0"），必须转数值后再比较
    out["max_phase"] = pd.to_numeric(out["max_phase"], errors="coerce")
    out.to_csv(TAB / "mi_gold_drugs.csv", index=False, encoding="utf-8-sig")

    # 无结构分子名单（结构缺失者排除出评估集，保留在名单内备查）
    no_struct = out[out["smiles"].isna() | (out["smiles"].astype(str).str.strip() == "")]
    no_struct[["molecule_chembl_id", "pref_name", "max_phase", "efo_terms"]].to_csv(
        TAB / "mi_gold_no_structure.csv", index=False, encoding="utf-8-sig")

    n_total = len(out)
    approved = int((out["max_phase"] == 4).sum())
    clinical = int(out["max_phase"].isin([1, 2, 3]).sum())
    other = n_total - approved - clinical
    n_smiles = int(out["smiles"].notna().sum())

    # efo 词条分布
    term_counts = Counter()
    for s in df["efo_term"].dropna():
        for t in str(s).split(";"):
            if t.strip():
                term_counts[t.strip()] += 1

    print(f"MI 药物 {n_total} 个（有 SMILES {n_smiles}）：approved {approved} / clinical {clinical} / "
          f"other {other}")

    log = ["# 10b MI 适应症药物（评估金标准）", "",
           f"- 生效查询方式：{used}",
           f"- 适应证记录数：{len(df)}，去重分子数：{n_total}",
           f"- 有 canonical_smiles：{n_smiles}；无结构：{len(no_struct)}"
           "（结构缺失者排除出评估集，保留在名单内备查）",
           f"- max_phase 分布：approved(4)={approved}，clinical(1-3)={clinical}，其他/未知={other}",
           "- 用途：仅作评估金标准，绝不进入训练集",
           "",
           "## 适应症词条（efo_term）分布", ""]
    for t, c in sorted(term_counts.items(), key=lambda x: -x[1]):
        log.append(f"- {t}: {c}")
    log.append("")
    (LOGS / "mi_drugs_log.md").write_text("\n".join(log), encoding="utf-8")
    print("[10b] 完成")


if __name__ == "__main__":
    main()