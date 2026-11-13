# -*- coding: utf-8 -*-
"""20 构建已上市药物候选池（ChEMBL max_phase=4 小分子）。

输入：ChEMBL molecule 接口（缓存 data/raw/chembl/approved/）
输出：results/tables/approved_pool.csv、logs/p8_log.md（漏斗部分）
运行：C:/Anaconda/Anaconda3/python.exe src/20_build_approved_pool.py

漏斗：原始数 -> 有结构 -> 标准化后 -> 去重后 -> 剔除（与训练正样本重叠 / 与 MI 金标准重叠 / 无结构）
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, TAB  # noqa: E402
from pool_utils import CACHE_DIR, descriptors, prepare_smiles, chembl_get  # noqa: E402

PAGE = 1000
CACHE_SUB = CACHE_DIR / "approved"


def main():
    CACHE_SUB.mkdir(parents=True, exist_ok=True)
    rows, offset, total = [], 0, None
    while True:
        data, _ = chembl_get("/molecule.json",
                             {"max_phase": 4, "molecule_type": "Small molecule",
                              "limit": PAGE, "offset": offset,
                              "only": "molecule_chembl_id,pref_name,molecule_structures,max_phase"},
                             cache_name=f"approved/approved_off{offset}.json")
        total = data["page_meta"]["total_count"]
        items = data.get("molecules") or []
        rows.extend(items)
        offset += PAGE
        if not items or offset >= total:
            break
    print(f"ChEMBL max_phase=4 小分子总数：{total}，拉取 {len(rows)}")

    recs, no_struct, unparsable = [], 0, 0
    for m in rows:
        st = (m.get("molecule_structures") or {}).get("canonical_smiles")
        if not st:
            no_struct += 1
            continue
        mol, reason = prepare_smiles(st)
        if mol is None:
            unparsable += 1
            continue
        d = descriptors(mol)
        d.update({"molecule_chembl_id": m["molecule_chembl_id"],
                  "pref_name": m.get("pref_name"),
                  "max_phase": m.get("max_phase")})
        recs.append(d)

    df = pd.DataFrame(recs)
    n_std = len(df)
    df = df.drop_duplicates("inchikey14", keep="first")
    n_dedup = len(df)
    print(f"有结构 {len(rows) - no_struct}，标准化后 {n_std}，InChIKey14 去重后 {n_dedup}")

    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    train_ik = set(act["inchikey14"])
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")
    gold_ik = set()
    for s in gold["smiles"].dropna():
        mol, _ = prepare_smiles(s)
        if mol is not None:
            gold_ik.add(descriptors(mol)["inchikey14"])

    n_train_overlap = int(df["inchikey14"].isin(train_ik).sum())
    n_gold_overlap = int(df["inchikey14"].isin(gold_ik).sum())
    keep = df[~df["inchikey14"].isin(train_ik | gold_ik)].copy()

    out = keep[["molecule_chembl_id", "pref_name", "smiles", "inchikey14", "mw", "clogp",
                "max_phase"]]
    out.to_csv(TAB / "approved_pool.csv", index=False, encoding="utf-8-sig")

    log = ["# P8 已上市药物重定位", "",
           "## 一、候选池漏斗", "",
           f"- ChEMBL max_phase=4 小分子总数：{total}",
           f"- 拉取记录：{len(rows)}",
           f"- 无结构：{no_struct}",
           f"- 标准化失败（无法解析/金属/重原子<8/MW>800）：{unparsable}",
           f"- 标准化后：{n_std}",
           f"- InChIKey14 去重后：{n_dedup}",
           f"- 与训练正样本（pool_actives）重叠剔除：{n_train_overlap}",
           f"- 与 MI 金标准重叠剔除：{n_gold_overlap}",
           f"- **最终候选池：{len(out)}**",
           ""]
    (LOGS / "p8_log.md").write_text("\n".join(log), encoding="utf-8")
    print(f"[20] 最终候选池 {len(out)}（训练重叠 {n_train_overlap}，MI 重叠 {n_gold_overlap}）")


if __name__ == "__main__":
    main()