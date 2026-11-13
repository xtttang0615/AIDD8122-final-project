# -*- coding: utf-8 -*-
"""15 Top 候选：每个主靶点取诱饵池模型概率 top-20，用 ChEMBL 补 pref_name。

输入：data/processed/ligand_features_ecfp4.npz、scaffold_splits.csv、pool_decoys.csv、
      results/tables/ligand_model_perf.csv（主模型名）
输出：results/tables/top_candidates.csv、logs/p7_experiments.md（追加）
运行：C:/Anaconda/Anaconda3/python.exe src/15_top_candidates.py

机制/适应症注释列留空，不凭记忆填写，由人工文献核对。
"""
from pathlib import Path
import sys
import datetime

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import LOGS, MAIN_TARGETS, TAB, load_features, load_splits, make_model  # noqa: E402
from pool_utils import chembl_get  # noqa: E402

TOP_N = 20


def main():
    ecfp, _ = load_features()
    splits = load_splits()
    perf = pd.read_csv(TAB / "ligand_model_perf.csv", encoding="utf-8-sig")
    main_model = perf.loc[perf["is_main"], "model"].iloc[0]
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")

    picks = []
    for target in MAIN_TARGETS:
        df = splits[splits["target"] == target].reset_index(drop=True)
        y = df["y"].to_numpy()
        X = np.vstack([ecfp[m] for m in df["molecule_chembl_id"]]).astype(np.float32)
        model = make_model(main_model, y)
        model.fit(X, y)

        d_ids = dec.loc[dec["target"] == target, "molecule_chembl_id"].tolist()
        Xd = np.vstack([ecfp[m] for m in d_ids]).astype(np.float32)
        p = model.predict_proba(Xd)[:, 1]
        order = np.argsort(-p)[:TOP_N]
        for rank, i in enumerate(order, 1):
            picks.append({"target": target, "rank": rank,
                          "molecule_chembl_id": d_ids[i], "prob": float(p[i])})

    ids = sorted({p["molecule_chembl_id"] for p in picks})
    names = {}
    for i in range(0, len(ids), 100):
        chunk = ids[i:i + 100]
        data, _ = chembl_get("/molecule.json",
                             {"molecule_chembl_id__in": ",".join(chunk), "limit": 100,
                              "only": "molecule_chembl_id,pref_name"})
        for m in data.get("molecules") or []:
            names[m["molecule_chembl_id"]] = m.get("pref_name")
    print(f"补名分子 {len(ids)} 个，命中 {sum(1 for v in names.values() if v)}")

    out = pd.DataFrame(picks)
    out["pref_name"] = out["molecule_chembl_id"].map(names)
    out["mechanism"] = ""            # 留空，人工文献核对
    out["indication_note"] = ""      # 留空，人工文献核对
    out = out[["target", "rank", "molecule_chembl_id", "pref_name", "prob",
               "mechanism", "indication_note"]]
    out.to_csv(TAB / "top_candidates.csv", index=False, encoding="utf-8-sig")

    with open(LOGS / "p7_experiments.md", "a", encoding="utf-8") as f:
        f.write(f"\n## {datetime.datetime.now().isoformat(timespec='seconds')} Top 候选（15）\n")
        f.write(f"- 主模型 {main_model}；每靶点取诱饵池概率 top-{TOP_N}；"
                f"ChEMBL 补名 {len(ids)} 个分子\n")
        f.write("- mechanism / indication_note 留空，由人工文献核对，未凭记忆填写\n\n")

    print(out.head(10).to_string(index=False))
    print("[15_top_candidates] 完成")


if __name__ == "__main__":
    main()