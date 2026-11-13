# -*- coding: utf-8 -*-
"""10d RDKit 标准化与合并：正样本 / 诱饵 / 清单 / 泄漏自检。

输入：data/interim/actives_raw.csv、data/interim/decoys_raw.csv、
      results/tables/mi_gold_drugs.csv、logs/actives_fetch_log.md、logs/decoy_log.md
输出：results/tables/pool_actives.csv、pool_decoys.csv、pool_manifest.csv、logs/pool_log.md
运行：C:/Anaconda/Anaconda3/python.exe src/10d_standardize.py

标准化：最大片段 -> 去盐 -> 电荷中和 -> canonical tautomer；
丢弃：无机/金属、重原子<8、MW>800、无法解析（逐类计数）；
去重：InChIKey14，保留 pChEMBL 最高者；单靶点正样本上限 800。
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
from pool_utils import INTERIM, LOGS, ROOT, descriptors, prepare_smiles  # noqa: E402

TAB = ROOT / "results" / "tables"
CAP_ACTIVES = 800
EXPLORATORY = {"PDK4"}  # 正样本数不足，骨架切分后测试集过小


def standardize_table(df, smiles_col="smiles"):
    kept, drop = [], Counter()
    for r in df.itertuples(index=False):
        mol, reason = prepare_smiles(getattr(r, smiles_col))
        if mol is None:
            drop[reason] += 1
            continue
        d = descriptors(mol)
        rec = {"target": r.target, "molecule_chembl_id": r.molecule_chembl_id}
        rec.update(d)
        if hasattr(r, "pchembl"):
            rec["pchembl"] = float(r.pchembl)
        kept.append(rec)
    return pd.DataFrame(kept), drop


def main():
    raw = pd.read_csv(INTERIM / "actives_raw.csv", encoding="utf-8-sig")
    dec_raw = pd.read_csv(INTERIM / "decoys_raw.csv", encoding="utf-8-sig")
    gold = pd.read_csv(TAB / "mi_gold_drugs.csv", encoding="utf-8-sig")

    act_std, act_drop = standardize_table(raw)
    dec_std, dec_drop = standardize_table(dec_raw)

    # 每条记录的原始活性数、空 SMILES 数（empty_smiles 单独计数，不与 parse_fail 混同）
    raw["is_empty"] = raw["smiles"].isna() | (raw["smiles"].astype(str).str.strip() == "")

    manifest, act_out, dec_out = [], [], []
    funnel = []

    for target, g in act_std.groupby("target"):
        n_rec = int((raw["target"] == target).sum())
        n_empty = int(raw.loc[raw["target"] == target, "is_empty"].sum())
        g = g.sort_values("pchembl", ascending=False).drop_duplicates("inchikey14", keep="first")
        n_dedup = len(g)
        capped = n_dedup > CAP_ACTIVES
        if capped:
            g = g.nlargest(CAP_ACTIVES, "pchembl")
        act_out.append(g)
        funnel.append({"target": target, "n_records": n_rec, "n_empty": n_empty,
                       "n_std": int((act_std["target"] == target).sum()),
                       "n_dedup": n_dedup, "n_final": len(g), "capped": capped})

    actives = pd.concat(act_out, ignore_index=True)
    actives = actives[["target", "molecule_chembl_id", "smiles", "inchikey14", "pchembl", "mw", "clogp"]]
    actives.to_csv(TAB / "pool_actives.csv", index=False, encoding="utf-8-sig")

    decoys = dec_std[["target", "molecule_chembl_id", "smiles", "inchikey14", "mw", "clogp"]]
    decoys.to_csv(TAB / "pool_decoys.csv", index=False, encoding="utf-8-sig")

    # 泄漏自检：MI 金标准与正样本的 InChIKey14 交集
    gold_keys = {}
    gold_drop = Counter()
    for r in gold.itertuples(index=False):
        mol, reason = prepare_smiles(r.smiles)
        if mol is None:
            gold_drop[reason] += 1
            continue
        gold_keys[r.molecule_chembl_id] = descriptors(mol)["inchikey14"]
    gold_set = set(gold_keys.values())

    leak_rows = []
    for target, g in actives.groupby("target"):
        inter = set(g["inchikey14"]) & gold_set
        leak_rows.append({"target": target, "n_actives": len(g), "n_overlap_gold": len(inter),
                          "overlap_inchikey14": ";".join(sorted(inter))})

    for f in funnel:
        t = f["target"]
        manifest.append({
            "target": t,
            "exploratory": t in EXPLORATORY,
            "n_actives_records": f["n_records"],
            "n_empty_smiles": f["n_empty"],
            "n_actives_standardized": f["n_std"],
            "n_actives_dedup": f["n_dedup"],
            "n_actives_final": f["n_final"],
            "capped_at_800": f["capped"],
            "n_decoys": int((decoys["target"] == t).sum()),
            "decoy_ratio": round(int((decoys["target"] == t).sum()) / max(f["n_final"], 1), 2),
        })
    man = pd.DataFrame(manifest).sort_values("n_actives_final", ascending=False)
    man.to_csv(TAB / "pool_manifest.csv", index=False, encoding="utf-8-sig")

    total_overlap = int(sum(r["n_overlap_gold"] for r in leak_rows))

    # ---------- 日志
    no_struct = gold[gold["smiles"].isna() | (gold["smiles"].astype(str).str.strip() == "")]
    term_counts = Counter()
    for s in gold["efo_terms"].dropna():
        for t in str(s).split(";"):
            if t.strip():
                term_counts[t.strip()] += 1

    log = ["# P6 分子池构建日志", "",
           "## 一、正样本漏斗（每靶点）", "",
           "| target | exploratory | 原始分子 | 空 SMILES | 标准化后 | InChIKey14 去重 | 最终正样本 | capped |",
           "|---|---|---|---|---|---|---|---|"]
    for f in funnel:
        log.append(f"| {f['target']} | {f['target'] in EXPLORATORY} | {f['n_records']} | {f['n_empty']} | "
                   f"{f['n_std']} | {f['n_dedup']} | {f['n_final']} | "
                   f"{'是（capped at 800）' if f['capped'] else '否'} |")
    capped_detail = [f for f in funnel if f["capped"]]
    if capped_detail:
        log.append("")
        for f in capped_detail:
            log.append(f"- cap 明细（{f['target']}）：{f['n_records']} → InChIKey14 去重 {f['n_dedup']} "
                       f"→ 取 pChEMBL 前 {CAP_ACTIVES} → {f['n_final']}")
    log += ["",
            f"- exploratory 靶点：{', '.join(sorted(EXPLORATORY))}"
            "（理由：正样本数不足，骨架切分后测试集过小）",
            "",
            "## 二、丢弃原因计数", "",
            f"- 正样本：{dict(act_drop) if act_drop else '无'}",
            f"- 诱饵：{dict(dec_drop) if dec_drop else '无'}",
            f"- MI 金标准：{dict(gold_drop) if gold_drop else '无'}",
            "- empty_smiles（接口未返回结构）与 parse_fail（RDKit 无法解析）分开计数",
            "",
            "## 三、诱饵", "",
            "| target | 正样本 | 诱饵数 | 比例 |",
            "|---|---|---|---|"]
    for r in man.itertuples(index=False):
        log.append(f"| {r.target} | {r.n_actives_final} | {r.n_decoys} | {r.decoy_ratio}:1 |")
    log += ["",
            "## 四、MI 药物金标准", "",
            f"- 漏斗：{len(gold)} 个药物 → 有结构 {int(gold['smiles'].notna().sum())} → "
            f"标准化后可用 {len(gold_keys)}",
            f"- 无结构 {len(no_struct)}（结构缺失者排除出评估集，保留在名单内备查）",
            f"- 查询方式：{gold['query_method'].iloc[0]}",
            f"- max_phase：approved(4)={int((gold['max_phase'] == 4).sum())}，"
            f"clinical(1-3)={int(gold['max_phase'].isin([1, 2, 3]).sum())}",
            "",
            "### efo_term 分布（按药物数）", ""]
    for t, c in sorted(term_counts.items(), key=lambda x: -x[1]):
        log.append(f"- {t}: {c}")
    log += ["",
            "## 五、泄漏自检（评估集 vs 训练正样本，InChIKey14）", ""]
    for r in leak_rows:
        log.append(f"- {r['target']}: 正样本 {r['n_actives']}，与 MI 金标准重叠 {r['n_overlap_gold']}")
    log += [f"- **评估集与训练集重叠数 = {total_overlap}，P7 将剔除**", ""]
    (LOGS / "pool_log.md").write_text("\n".join(log), encoding="utf-8")

    print(f"正样本 {len(actives)}，诱饵 {len(decoys)}，MI 药物 {len(gold)}")
    print(f"泄漏重叠总数 = {total_overlap}")
    print(man.to_string(index=False))
    print("[10d] 完成")


if __name__ == "__main__":
    main()