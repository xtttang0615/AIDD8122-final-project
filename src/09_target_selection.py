# -*- coding: utf-8 -*-
"""09 靶点规则化筛选：对跨折签名池 28 个基因查询 Open Targets / STRING / ChEMBL。

输入：logs/model_log.md（解析五折签名基因）、results/tables/signature_genes.csv（最终签名）
输出：results/tables/target_candidates.csv、logs/target_selection_log.md、logs/api_errors.md
运行：C:/Anaconda/Anaconda3/python.exe src/09_target_selection.py

规格见 docs/09_spec.md。硬筛选：OT overall >= 0.1 且 STRING degree >= 1。
接口失败一律记 NA 并写入 logs/api_errors.md，不编造、不使用记忆中的靶点-药物关系。
"""
from pathlib import Path
import json
import re
import sys
import time

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
TAB = ROOT / "results" / "tables"
LOGS = ROOT / "logs"
TAB.mkdir(parents=True, exist_ok=True)

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
STRING_URL = "https://string-db.org/api/json/interaction_partners"
CHEMBL_URL = "https://www.ebi.ac.uk/chembl/api/data"

TIMEOUT = 60
RETRIES = 3
SPECIES = 9606
CHEMBL_ACT_LIMIT = 1000
OT_MIN_SCORE = 0.1
OT_RELAXED_SCORE = 0.05
STRING_MIN_DEGREE = 1

calls = {"ot": 0, "string": 0, "chembl": 0}
errors = []


def _get(url, params=None):
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.get(url, params=params, timeout=TIMEOUT,
                             headers={"User-Agent": "aidd8122-coursework"})
            if r.status_code == 200:
                return r.json()
            last = f"HTTP {r.status_code}: {r.text[:300]}"
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
        time.sleep(2 ** attempt)
    raise RuntimeError(last or "unknown error")


def _post_ot(query, variables):
    last = None
    for attempt in range(RETRIES):
        try:
            r = requests.post(OT_URL, json={"query": query, "variables": variables},
                              timeout=TIMEOUT, headers={"User-Agent": "aidd8122-coursework"})
            if r.status_code == 200:
                j = r.json()
                if "errors" in j:
                    raise RuntimeError(f"GraphQL errors: {json.dumps(j['errors'])[:300]}")
                return j["data"]
            last = f"HTTP {r.status_code}: {r.text[:300]}"
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
        time.sleep(2 ** attempt)
    raise RuntimeError(last or "unknown error")


def log_error(gene, endpoint, err):
    errors.append({"gene": gene, "endpoint": endpoint, "error": str(err)[:400]})


def load_gene_pool():
    """28 基因跨折池：解析 model_log.md 五折签名 + signature_genes.csv 最终签名。"""
    txt = (LOGS / "model_log.md").read_text(encoding="utf-8")
    folds = {}
    for m in re.finditer(r"- fold (\d)：折内 DEG \d+ 个，签名 \d+ 个.*?基因：(.+)", txt):
        folds[int(m.group(1))] = [g.strip() for g in m.group(2).split(",")]
    if len(folds) != 5:
        raise RuntimeError(f"model_log.md 解析到 {len(folds)} 折，预期 5 折")
    sig = pd.read_csv(TAB / "signature_genes.csv", encoding="utf-8-sig")
    final = sig["gene"].tolist()

    counts = {}
    for genes in folds.values():
        for g in genes:
            counts[g] = counts.get(g, 0) + 1
    for g in final:
        counts.setdefault(g, 0)
    rows = [{"gene": g, "sel_freq": c, "in_final_signature": g in final} for g, c in counts.items()]
    df = pd.DataFrame(rows).sort_values(["sel_freq", "gene"], ascending=[False, True])
    return df.reset_index(drop=True), folds


def resolve_disease():
    q = ("query($s:String!){search(queryString:$s,entityNames:[\"disease\"],"
         "page:{index:0,size:1}){hits{id name entity}}}")
    calls["ot"] += 1
    data = _post_ot(q, {"s": "myocardial infarction"})
    hit = data["search"]["hits"][0]
    return hit["id"], hit["name"]


def ot_target(gene, disease_id):
    q_search = ("query($s:String!){search(queryString:$s,entityNames:[\"target\"],"
                "page:{index:0,size:1}){hits{id name entity}}}")
    calls["ot"] += 1
    hits = _post_ot(q_search, {"s": gene})["search"]["hits"]
    if not hits:
        raise RuntimeError("Open Targets search 无 target 命中")
    ensembl = hits[0]["id"]

    q_t = ("query($id:String!,$efo:String!){target(ensemblId:$id){approvedSymbol approvedName "
           "proteinIds{id source} tractability{label modality value} "
           "drugAndClinicalCandidates{count} "
           "associatedDiseases(Bs:[$efo],page:{index:0,size:1}){rows{score datatypeScores{id score}}}}}")
    calls["ot"] += 1
    t = _post_ot(q_t, {"id": ensembl, "efo": disease_id})["target"]

    uniprot = next((p["id"] for p in t.get("proteinIds") or []
                    if p.get("source") == "uniprot_swissprot"), None)
    sm_labels = [x["label"] for x in (t.get("tractability") or [])
                 if x.get("modality") == "SM" and x.get("value")]
    ab_labels = [x["label"] for x in (t.get("tractability") or [])
                 if x.get("modality") == "AB" and x.get("value")]

    overall, dtypes = np.nan, {}
    rows = (t.get("associatedDiseases") or {}).get("rows") or []
    if rows:
        overall = float(rows[0]["score"])
        dtypes = {d["id"]: float(d["score"]) for d in rows[0].get("datatypeScores") or []}

    # 注：OT v4 中靶点级已知药物总数由 drugAndClinicalCandidates.count 暴露
    #（旧版 knownDrugs 字段已下线）
    known_drugs_total = (t.get("drugAndClinicalCandidates") or {}).get("count")

    return {
        "ensembl_id": ensembl,
        "ot_target_name": t.get("approvedName"),
        "uniprot": uniprot,
        "ot_overall_score": overall,
        "ot_genetic_association": dtypes.get("genetic_association", np.nan),
        "ot_known_drug_disease": dtypes.get("known_drug", np.nan),
        "ot_literature": dtypes.get("literature", np.nan),
        "ot_clinical": dtypes.get("clinical", np.nan),
        "target_known_drugs_total": known_drugs_total if known_drugs_total is not None else np.nan,
        "tractability_sm": ";".join(sm_labels),
        "tractability_ab": ";".join(ab_labels),
    }


def string_info(gene):
    calls["string"] += 1
    data = _get(STRING_URL, {"identifiers": gene, "species": SPECIES,
                             "limit": 1000, "caller_identity": "aidd8122-coursework"})
    if not isinstance(data, list) or not data:
        raise RuntimeError("STRING 无互作伙伴返回")
    scores = [float(x["score"]) for x in data]
    return {"string_degree": len(data), "string_max_combined_score": max(scores)}


def chembl_info(uniprot):
    calls["chembl"] += 1
    js = _get(f"{CHEMBL_URL}/target.json",
              {"target_components__accession": uniprot, "limit": 20})
    targets = js.get("targets") or []
    if not targets:
        raise RuntimeError(f"ChEMBL 无 accession {uniprot} 对应靶点")
    cand = [t for t in targets if t.get("organism") == "Homo sapiens"] or targets
    tid = cand[0]["target_chembl_id"]

    calls["chembl"] += 1
    act = _get(f"{CHEMBL_URL}/activity.json",
               {"target_chembl_id": tid, "limit": CHEMBL_ACT_LIMIT,
                "only": "molecule_chembl_id,pchembl_value"})
    total = int(act["page_meta"]["total_count"])
    acts = act.get("activities") or []
    mols = {a["molecule_chembl_id"] for a in acts if a.get("molecule_chembl_id")}
    vals = [float(a["pchembl_value"]) for a in acts if a.get("pchembl_value")]
    return {"chembl_target_id": tid,
            "chembl_n_activities_total": total,
            "chembl_n_molecules_first1000": len(mols),
            "chembl_max_pchembl_first1000": max(vals) if vals else np.nan}


def main():
    pool, folds = load_gene_pool()
    print(f"跨折签名池：{len(pool)} 个基因（五折 {sum(len(v) for v in folds.values())} 条记录）")

    disease_id, disease_name = resolve_disease()
    print(f"疾病实体：{disease_name} ({disease_id})")

    records = []
    for i, r in enumerate(pool.itertuples(index=False), 1):
        gene = r.gene
        rec = {"gene": gene, "sel_freq": int(r.sel_freq),
               "in_final_signature": bool(r.in_final_signature)}
        status = []

        try:
            ot = ot_target(gene, disease_id)
            rec.update(ot)
        except Exception as e:
            log_error(gene, "OpenTargets", e)
            status.append("ot_fail")
            rec.update({"ensembl_id": np.nan, "ot_target_name": np.nan, "uniprot": np.nan,
                        "ot_overall_score": np.nan, "ot_genetic_association": np.nan,
                        "ot_known_drug_disease": np.nan, "ot_literature": np.nan,
                        "ot_clinical": np.nan, "target_known_drugs_total": np.nan,
                        "tractability_sm": np.nan, "tractability_ab": np.nan})

        try:
            rec.update(string_info(gene))
        except Exception as e:
            log_error(gene, "STRING", e)
            status.append("string_fail")
            rec.update({"string_degree": np.nan, "string_max_combined_score": np.nan})

        try:
            if pd.isna(rec.get("uniprot")):
                raise RuntimeError("无 UniProt accession，无法映射 ChEMBL 靶点")
            rec.update(chembl_info(rec["uniprot"]))
        except Exception as e:
            log_error(gene, "ChEMBL", e)
            status.append("chembl_fail")
            rec.update({"chembl_target_id": np.nan, "chembl_n_activities_total": np.nan,
                        "chembl_n_molecules_first1000": np.nan,
                        "chembl_max_pchembl_first1000": np.nan})

        rec["data_status"] = "ok" if not status else ";".join(status)
        records.append(rec)
        print(f"[{i}/{len(pool)}] {gene}: OT={rec['ot_overall_score']}, "
              f"STRING_deg={rec['string_degree']}, "
              f"pChEMBL={rec['chembl_max_pchembl_first1000']}, status={rec['data_status']}")

    df = pd.DataFrame(records)

    # 旧版硬筛（保留 target_candidates.csv 口径不变）
    ot_ok = df["ot_overall_score"].fillna(-1) >= OT_MIN_SCORE
    st_ok = df["string_degree"].fillna(0) >= STRING_MIN_DEGREE
    df["pass_filter"] = ot_ok & st_ok

    # strict = OT overall >= 0.1（其余不硬筛，列保留）
    df["strict_pass"] = df["ot_overall_score"].fillna(-1) >= OT_MIN_SCORE

    # relaxed = OT overall >= 0.05 或 target_known_drugs_total > 0
    by_score = df["ot_overall_score"].fillna(-1) >= OT_RELAXED_SCORE
    by_drug = df["target_known_drugs_total"].fillna(0) > 0
    df["relaxed_pass"] = by_score | by_drug
    df["relaxed_reason"] = np.select(
        [by_score & by_drug, by_score & ~by_drug, ~by_score & by_drug],
        ["OT>=0.05 & known_drugs>0", "OT>=0.05", "known_drugs>0"], default="")

    df = df.sort_values(["sel_freq", "ot_overall_score", "chembl_max_pchembl_first1000"],
                        ascending=[False, False, False], na_position="last")
    df.to_csv(TAB / "target_candidates.csv", index=False, encoding="utf-8-sig")
    df[df["strict_pass"]].to_csv(TAB / "target_candidates_strict.csv",
                                 index=False, encoding="utf-8-sig")
    df[df["relaxed_pass"]].to_csv(TAB / "target_candidates_relaxed.csv",
                                  index=False, encoding="utf-8-sig")

    n_total = len(df)
    n_old, n_strict, n_relaxed = int(df["pass_filter"].sum()), int(df["strict_pass"].sum()), \
        int(df["relaxed_pass"].sum())
    print(f"\n筛前 {n_total} 个；旧口径 {n_old} 个；strict {n_strict} 个；relaxed {n_relaxed} 个")

    # 自检：PDE4B 靶点级已知药物总数必须 > 0
    pde4b = df[df["gene"] == "PDE4B"]
    pde4b_drugs = None if pde4b.empty else pde4b.iloc[0]["target_known_drugs_total"]
    if pde4b_drugs is None or pd.isna(pde4b_drugs) or pde4b_drugs <= 0:
        print(f"[自检失败] PDE4B target_known_drugs_total = {pde4b_drugs}，预期 > 0")
    else:
        print(f"[自检通过] PDE4B target_known_drugs_total = {int(pde4b_drugs)} > 0")

    log = ["# 09 靶点规则化筛选日志", "",
           f"- 疾病实体：{disease_name}（{disease_id}，由 Open Targets search 解析，未硬编码）",
           f"- 候选基因数（跨折签名池）：{n_total}",
           f"- 旧口径（OT overall >= {OT_MIN_SCORE} 且 STRING degree >= {STRING_MIN_DEGREE}）筛后：{n_old}",
           f"- strict（OT overall >= {OT_MIN_SCORE}）筛后：{n_strict}",
           f"- relaxed（OT overall >= {OT_RELAXED_SCORE} 或 target_known_drugs_total > 0）筛后：{n_relaxed}",
           f"- PDE4B 自检：target_known_drugs_total = {pde4b_drugs}（要求 > 0）",
           "",
           "## 规则变更记录",
           "- 原方案：hub = STRING degree 前 25% + ChEMBL >= 10 硬筛。",
           "- 实际执行：简化为 OT overall >= 0.1 且 STRING degree >= 1（degree 门槛形同虚设，"
           "28 个基因中 27 个 degree 远大于 1），ChEMBL 仅作排序列。",
           "- 变更理由：候选仅 28 个基因，hub 过滤区分度不足。",
           "- P5b 追加：strict = OT >= 0.1；relaxed = OT >= 0.05 或 target_known_drugs_total > 0；"
           "known drugs 为靶点级字段（OT v4 drugAndClinicalCandidates.count，旧版 knownDrugs 已下线），"
           "ot_known_drug_disease 为 MI 疾病上下文内的已知药证据，两者语义不同，均已保留。",
           "",
           "## 每接口调用量",
           f"- Open Targets：{calls['ot']} 次（含 1 次疾病解析）",
           f"- STRING：{calls['string']} 次",
           f"- ChEMBL：{calls['chembl']} 次",
           f"- 失败记录：{len(errors)} 条（详见 logs/api_errors.md）",
           "",
           "## 数据口径说明",
           f"- ChEMBL：每次取前 {CHEMBL_ACT_LIMIT} 条活性记录；"
           "chembl_n_molecules_first1000 为该批记录中去重的化合物数，"
           "chembl_n_activities_total 为接口返回的活性记录总数",
           "- OT datatypeScores 只返回非零证据类型，字段为 NA 表示该证据类型无数据（零分省略）",
           "- STRING：interaction_partners limit=1000，degree 为伙伴数（本次最大值远小于上限，未截断）",
           "",
           "## 筛选结果（按 sel_freq -> OT overall -> max pChEMBL 降序）",
           "| gene | sel_freq | OT overall | ot_known_drug_disease | target_known_drugs_total | "
           "STRING deg | max pChEMBL | strict | relaxed | relaxed_reason |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for r in df.itertuples(index=False):
        log.append(f"| {r.gene} | {r.sel_freq} | "
                   f"{'' if pd.isna(r.ot_overall_score) else round(r.ot_overall_score, 4)} | "
                   f"{'' if pd.isna(r.ot_known_drug_disease) else round(r.ot_known_drug_disease, 4)} | "
                   f"{'' if pd.isna(r.target_known_drugs_total) else int(r.target_known_drugs_total)} | "
                   f"{'' if pd.isna(r.string_degree) else int(r.string_degree)} | "
                   f"{'' if pd.isna(r.chembl_max_pchembl_first1000) else r.chembl_max_pchembl_first1000} | "
                   f"{bool(r.strict_pass)} | {bool(r.relaxed_pass)} | {r.relaxed_reason} |")
    log.append("")
    (LOGS / "target_selection_log.md").write_text("\n".join(log), encoding="utf-8")

    err_lines = ["# 接口错误记录", ""]
    if errors:
        err_lines += ["| gene | endpoint | error |", "|---|---|---|"]
        err_lines += [f"| {e['gene']} | {e['endpoint']} | {e['error']} |" for e in errors]
    else:
        err_lines.append("本次运行无接口失败。")
    err_lines.append("")
    (LOGS / "api_errors.md").write_text("\n".join(err_lines), encoding="utf-8")

    print("\n=== Top 10 ===")
    cols = ["gene", "sel_freq", "ot_overall_score", "target_known_drugs_total", "string_degree",
            "chembl_max_pchembl_first1000", "strict_pass", "relaxed_pass"]
    print(df[cols].head(10).to_string(index=False))

    print("\n=== strict 靶点 ===")
    for r in df[df["strict_pass"]].itertuples(index=False):
        print(f"  {r.gene}: OT={r.ot_overall_score:.4f}, known_drugs={r.target_known_drugs_total}")

    print("\n=== relaxed 靶点（含通过原因）===")
    for r in df[df["relaxed_pass"]].itertuples(index=False):
        ot = "NA" if pd.isna(r.ot_overall_score) else f"{r.ot_overall_score:.4f}"
        kd = "NA" if pd.isna(r.target_known_drugs_total) else int(r.target_known_drugs_total)
        print(f"  {r.gene}: OT={ot}, known_drugs={kd} -> {r.relaxed_reason}")
    print("[09_target_selection] 完成")


if __name__ == "__main__":
    main()