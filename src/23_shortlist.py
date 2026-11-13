# -*- coding: utf-8 -*-
"""23 重定位短名单：按 prob_mean 取 top 50（排除 MI 金标准），另出各靶点 top 15。

输入：results/tables/approved_scores.csv、approved_pool.csv（取 max_phase）
输出：results/tables/repurposing_shortlist.csv、shortlist_by_target.csv、logs/p8_log.md（追加）
运行：C:/Anaconda/Anaconda3/python.exe src/23_shortlist.py

机制/适应症/文献列留空，禁止凭记忆填写。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, MAIN_TARGETS, TAB  # noqa: E402


def train_keys(target):
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    return (set(act.loc[act["target"] == target, "inchikey14"])
            | set(dec.loc[dec["target"] == target, "inchikey14"]))


def main():
    sc = pd.read_csv(TAB / "approved_scores.csv", encoding="utf-8-sig")
    pool = pd.read_csv(TAB / "approved_pool.csv", encoding="utf-8-sig")

    # 聚合排名要求分子对三个模型均未见过（训练正例与负例都算）
    seen = set()
    for t in MAIN_TARGETS:
        seen |= train_keys(t)
    n_before = len(sc[sc["is_mi_gold"] == 0])
    cand = sc[(sc["is_mi_gold"] == 0) & (~sc["inchikey14"].isin(seen))].copy()
    print(f"短名单候选：{n_before} -> {len(cand)}（剔除与任一模型训练数据重叠 {n_before - len(cand)}）")
    cand = cand.merge(pool[["molecule_chembl_id", "max_phase"]],
                      on="molecule_chembl_id", how="left")
    cand = cand.sort_values("prob_mean", ascending=False)

    top = cand.head(50).reset_index(drop=True)
    out = pd.DataFrame({
        "rank": range(1, len(top) + 1),
        "pref_name": top["pref_name"],
        "chembl_id": top["molecule_chembl_id"],
        "max_phase": top["max_phase"],
        "prob_PDE4B": top["prob_PDE4B"],
        "prob_PFKFB3": top["prob_PFKFB3"],
        "prob_EPAS1": top["prob_EPAS1"],
        "prob_mean": top["prob_mean"],
        "mechanism": "",
        "indication": "",
        "literature": "",
    })
    out.to_csv(TAB / "repurposing_shortlist.csv", index=False, encoding="utf-8-sig")

    by_t = []
    for t in MAIN_TARGETS:
        sub = cand.sort_values(f"prob_{t}", ascending=False).head(15)
        for i, r in enumerate(sub.itertuples(index=False), 1):
            by_t.append({"target": t, "rank": i, "pref_name": r.pref_name,
                         "chembl_id": r.molecule_chembl_id, "max_phase": r.max_phase,
                         "prob": getattr(r, f"prob_{t}")})
    pd.DataFrame(by_t).to_csv(TAB / "shortlist_by_target.csv", index=False, encoding="utf-8-sig")

    with open(LOGS / "p8_log.md", "a", encoding="utf-8") as f:
        f.write("## 四、重定位短名单（top 20 预览，已排除 MI 金标准）\n\n")
        f.write("| rank | pref_name | chembl_id | prob_mean |\n|---|---|---|---|\n")
        for r in out.head(20).itertuples(index=False):
            f.write(f"| {r.rank} | {r.pref_name} | {r.chembl_id} | {r.prob_mean:.4f} |\n")
        f.write("\n机制/适应症/文献列一律留空，待人工核对。\n\n")
        f.write("> 模型饱和，排序仅用于相对优先级，不作活性预测。\n")

    print(out.head(20)[["rank", "pref_name", "chembl_id", "prob_mean"]].to_string(index=False))
    print("[23_shortlist] 完成")


if __name__ == "__main__":
    main()