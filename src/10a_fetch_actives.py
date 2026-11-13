# -*- coding: utf-8 -*-
"""10a 拉取 relaxed 靶点的 ChEMBL 活性分子（正样本原始数据）。

输入：results/tables/target_candidates_relaxed.csv（7 靶点，含 chembl_target_id）
输出：data/interim/actives_raw.csv、logs/actives_fetch_log.md、logs/api_requests.log
运行：C:/Anaconda/Anaconda3/python.exe src/10a_fetch_actives.py

筛选：pchembl_value >= 6 且 standard_type in (IC50,Ki,Kd,EC50)；同分子多条取 max pChEMBL。
无配体靶点（chembl_target_id 为空）标记 no_ligand_data，不进建模。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pool_utils import (INTERIM, LOGS, ROOT, chembl_get, log_request_note)  # noqa: E402

TAB = ROOT / "results" / "tables"
PAGE = 1000
STD_TYPES = "IC50,Ki,Kd,EC50"
PCH_MAX_CHUNK = 100


def fetch_structures(ids):
    """批量取 canonical_smiles（每批 100 个）。"""
    out = {}
    ids = sorted(set(ids))
    for i in range(0, len(ids), PCH_MAX_CHUNK):
        chunk = ids[i:i + PCH_MAX_CHUNK]
        data, _ = chembl_get("/molecule.json",
                             {"molecule_chembl_id__in": ",".join(chunk), "limit": PCH_MAX_CHUNK,
                              "only": "molecule_chembl_id,molecule_structures"})
        for m in data.get("molecules") or []:
            st = m.get("molecule_structures") or {}
            out[m["molecule_chembl_id"]] = st.get("canonical_smiles")
    return out


def fetch_target_actives(tid):
    agg, n_records, offset = {}, 0, 0
    while True:
        data, _ = chembl_get("/activity.json",
                             {"target_chembl_id": tid, "pchembl_value__gte": 6,
                              "standard_type__in": STD_TYPES, "limit": PAGE, "offset": offset,
                              "only": "molecule_chembl_id,pchembl_value"},
                             cache_name=f"act_{tid}_off{offset}.json")
        total = data["page_meta"]["total_count"]
        acts = data.get("activities") or []
        n_records += len(acts)
        for a in acts:
            mid, pv = a.get("molecule_chembl_id"), a.get("pchembl_value")
            if mid and pv is not None:
                agg[mid] = max(agg.get(mid, 0.0), float(pv))
        offset += PAGE
        if not acts or offset >= total:
            break
    return agg, total, n_records


def main():
    cand = pd.read_csv(TAB / "target_candidates_relaxed.csv", encoding="utf-8-sig")
    log_request_note("10a 开始拉取活性数据")

    rows, logs = [], ["# 10a 活性数据拉取日志", ""]
    no_ligand = []

    for r in cand.itertuples(index=False):
        tid = r.chembl_target_id
        if not isinstance(tid, str) or not tid.strip() or tid == "nan":
            no_ligand.append(r.gene)
            logs.append(f"- {r.gene}: no_ligand_data（无 ChEMBL target id）")
            continue

        agg, total, n_records = fetch_target_actives(tid)
        if not agg:
            no_ligand.append(f"{r.gene}（有 target id {tid}，但 0 条合格活性）")
            logs.append(f"- {r.gene} ({tid}): no_ligand_data（0 条满足 pChEMBL>=6 的活性）")
            print(f"{r.gene} ({tid}): 0 条合格活性 -> no_ligand_data")
            continue
        structs = fetch_structures(agg.keys())
        n_struct = sum(1 for v in structs.values() if v)
        for mid, pv in agg.items():
            rows.append({"target": r.gene, "chembl_target_id": tid, "molecule_chembl_id": mid,
                         "smiles": structs.get(mid), "pchembl": pv})
        print(f"{r.gene} ({tid}): 活性记录 {n_records}/{total}，去重分子 {len(agg)}，有结构 {n_struct}")
        logs.append(f"- {r.gene} ({tid}): 过滤后活性记录 {n_records}（接口 total_count={total}），"
                    f"去重分子 {len(agg)}，其中有 canonical_smiles {n_struct}")

    df = pd.DataFrame(rows)
    df.to_csv(INTERIM / "actives_raw.csv", index=False, encoding="utf-8-sig")

    logs += ["", "## 无配体靶点（no_ligand_data，不进建模）",
             "- " + (", ".join(no_ligand) if no_ligand else "无"), ""]
    (LOGS / "actives_fetch_log.md").write_text("\n".join(logs), encoding="utf-8")

    print(f"\n[10a] 共 {len(df)} 条正样本记录，覆盖 {df['target'].nunique()} 个靶点；"
          f"no_ligand_data: {no_ligand}")


if __name__ == "__main__":
    main()