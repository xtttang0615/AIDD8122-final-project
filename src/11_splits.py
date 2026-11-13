# -*- coding: utf-8 -*-
"""11 骨架切分：Bemis-Murcko 骨架分组，生成并固化切分表。

输入：results/tables/pool_actives.csv、pool_decoys.csv
输出：data/processed/scaffold_splits.csv、logs/split_log.md
运行：C:/Anaconda/Anaconda3/python.exe src/11_splits.py

主切分：StratifiedGroupKFold(5, shuffle, random_state=42)，group = Murcko 骨架 canonical SMILES；
无骨架（全链分子）用 inchikey14 代替。PDK4：单次 80/20 骨架切分（seed=42），标 exploratory。
本表一次生成后固化，12/13 全部复用，不得重新生成。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import ALL_TARGETS, LOGS, MAIN_TARGETS, PROC, SPLIT_CSV, TAB  # noqa: E402

RDLogger.DisableLog("rdApp.*")
SEED = 42


def scaffold_of(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return ""
    try:
        return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    except Exception:
        return ""


def main():
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")

    out, logs = [], ["# P7 骨架切分日志", "",
                     "切分键：RDKit Bemis-Murcko 骨架 canonical SMILES；"
                     "无骨架分子用 InChIKey14 代替（记为 no_scaffold）",
                     f"主切分：StratifiedGroupKFold(5, shuffle=True, random_state={SEED})",
                     "PDK4：单次 80/20 骨架切分（seed=42），exploratory=True",
                     "本切分表已固化，后续模型与消融全部复用，不得重新生成", ""]

    for target in ALL_TARGETS:
        a = act[act["target"] == target][["molecule_chembl_id", "smiles", "inchikey14"]].copy()
        a["y"] = 1
        a["target"] = target
        d = dec[dec["target"] == target][["molecule_chembl_id", "smiles", "inchikey14"]].copy()
        d["y"] = 0
        d["target"] = target
        df = pd.concat([a, d], ignore_index=True)
        df["scaffold"] = df["smiles"].map(scaffold_of)
        no_scaf = int((df["scaffold"] == "").sum())
        df["group_key"] = np.where(df["scaffold"] == "", "ICK_" + df["inchikey14"], df["scaffold"])

        if target == "PDK4":
            rng = np.random.default_rng(SEED)
            groups = df["group_key"].unique()
            rng.shuffle(groups)
            n_test = max(1, int(round(len(groups) * 0.2)))
            test_set = set(groups[:n_test])
            df["fold"] = np.where(df["group_key"].isin(test_set), 1, 0)
            scheme = "holdout80_20"
        else:
            df["fold"] = -1
            sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
            for k, (_, te) in enumerate(sgkf.split(df, df["y"], groups=df["group_key"])):
                df.loc[te, "fold"] = k
            scheme = "cv5_stratified_group"

        df["scheme"] = scheme
        out.append(df[["target", "molecule_chembl_id", "y", "scaffold", "group_key",
                       "fold", "scheme"]])

        logs.append(f"## {target}"
                    + ("（exploratory）" if target == "PDK4" else ""))
        logs.append(f"- 样本：{len(df)}（正 {int(df.y.sum())} / 负 {int((1 - df.y).sum())}）")
        logs.append(f"- 骨架数：{df['scaffold'].nunique()}；无骨架分子：{no_scaf}；"
                    f"分组键数：{df['group_key'].nunique()}")
        if scheme == "holdout80_20":
            tr, te = df[df.fold == 0], df[df.fold == 1]
            logs.append(f"- 80/20：train {len(tr)}（正 {int(tr.y.sum())}）/ test {len(te)}（正 {int(te.y.sum())}）")
        else:
            for k in range(5):
                f = df[df.fold == k]
                logs.append(f"- fold {k}: n={len(f)}，正 {int(f.y.sum())}，负 {int((1 - f.y).sum())}")
        logs.append("")

    res = pd.concat(out, ignore_index=True)
    res["fold"] = res["fold"].astype(int)
    PROC.mkdir(parents=True, exist_ok=True)
    res.to_csv(SPLIT_CSV, index=False, encoding="utf-8-sig")
    (LOGS / "split_log.md").write_text("\n".join(logs), encoding="utf-8")

    print(res.groupby(["target", "fold"]).agg(n=("y", "size"), pos=("y", "sum")).to_string())
    print(f"[11_splits] 已固化 {SPLIT_CSV}")


if __name__ == "__main__":
    main()