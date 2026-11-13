# -*- coding: utf-8 -*-
"""10c 构建诱饵：从 ChEMBL Small molecule 抽样池按性质匹配，为每个靶点独立配 10:1 诱饵。

输入：data/interim/actives_raw.csv、ChEMBL molecule 接口（molecule_type=Small molecule）
输出：data/interim/decoys_raw.csv、logs/decoy_log.md
运行：C:/Anaconda/Anaconda3/python.exe src/10c_build_decoys.py

规则：seed=42；每个正样本配 10 个诱饵；诱饵与该靶点全部正样本 ECFP4 Tanimoto < 0.4；
性质窗口 MW±50、cLogP±1、HBD±1、HBA±1、RotB±2、Rings±1；
每靶点最多 20 轮，抽不满按实际数量记录，不硬凑；诱饵不跨靶点共享。
"""
from pathlib import Path
import sys

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

import numpy as np
import pandas as pd
from rdkit import DataStructs

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pool_utils import (INTERIM, LOGS, chembl_get, descriptors, log_request_note,  # noqa: E402
                        morgan_fp, prepare_smiles)

PAGE = 1000
N_PAGES = 30          # 30 x 1000 = 30000 个 Small molecule
SEED = 42
DECOYS_PER_ACTIVE = 10
MAX_ROUNDS = 20
TANIMOTO_MAX = 0.4
CAP_ACTIVES = 800
EXPLORATORY = {"PDK4"}  # 正样本数不足，骨架切分后测试集过小

WINDOW = {"mw": 50.0, "clogp": 1.0, "hbd": 1, "hba": 1, "rotb": 2, "rings": 1}


def build_actives():
    raw = pd.read_csv(INTERIM / "actives_raw.csv", encoding="utf-8-sig")
    recs = []
    for r in raw.itertuples(index=False):
        # 描述符与 Tanimoto 一律基于标准化后的正样本（parse -> standardize_mol -> 性质过滤）
        mol, _reason = prepare_smiles(r.smiles)
        if mol is None:
            continue
        d = descriptors(mol)
        d.update({"target": r.target, "molecule_chembl_id": r.molecule_chembl_id,
                  "pchembl": float(r.pchembl), "fp": morgan_fp(mol)})
        recs.append(d)
    df = pd.DataFrame(recs)
    # InChIKey14 去重，保留 pChEMBL 最高
    df = df.sort_values("pchembl", ascending=False).drop_duplicates("inchikey14", keep="first")
    capped = []
    parts = []
    for t, g in df.groupby("target"):
        if len(g) > CAP_ACTIVES:
            capped.append((t, len(g)))
            g = g.nlargest(CAP_ACTIVES, "pchembl")
        parts.append(g)
    df = pd.concat(parts, ignore_index=True)
    return df, capped


def build_decoy_pool():
    rows = []
    for page in range(N_PAGES):
        data, _ = chembl_get("/molecule.json",
                             {"molecule_type": "Small molecule", "limit": PAGE,
                              "offset": page * PAGE,
                              "only": "molecule_chembl_id,molecule_structures"},
                             cache_name=f"smallmol_off{page * PAGE}.json")
        for m in data.get("molecules") or []:
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if not smi:
                continue
            mol, _ = prepare_smiles(smi)
            if mol is None:
                continue
            d = descriptors(mol)
            d.update({"molecule_chembl_id": m["molecule_chembl_id"], "fp": morgan_fp(mol)})
            rows.append(d)
        print(f"  抽样池 page {page + 1}/{N_PAGES} 累计可用 {len(rows)}")
    df = pd.DataFrame(rows).drop_duplicates("inchikey14", keep="first").reset_index(drop=True)
    return df


def main():
    log_request_note("10c 开始构建诱饵")
    actives, capped = build_actives()
    print(f"标准化后正样本：{len(actives)} 条，覆盖 {actives['target'].nunique()} 个靶点")
    for t, n in capped:
        print(f"  [capped at {CAP_ACTIVES}] {t}: 去重后 {n} -> {CAP_ACTIVES}")

    print("拉取 Small molecule 抽样池 ...")
    pool = build_decoy_pool()
    print(f"抽样池可用分子：{len(pool)}")

    pool_fps = pool["fp"].tolist()
    p_mw = pool["mw"].to_numpy()
    p_logp = pool["clogp"].to_numpy()
    p_hbd = pool["hbd"].to_numpy()
    p_hba = pool["hba"].to_numpy()
    p_rotb = pool["rotb"].to_numpy()
    p_rings = pool["rings"].to_numpy()

    out_rows, logs = [], ["# 10c 诱饵构建日志", "",
                          f"- 抽样池来源：ChEMBL molecule_type=Small molecule，前 {N_PAGES * PAGE} 条"
                          "（按 offset 顺序取，非全库随机抽样）",
                          f"- seed={SEED}，每正样本 {DECOYS_PER_ACTIVE} 个诱饵，"
                          f"Tanimoto < {TANIMOTO_MAX}，每靶点最多 {MAX_ROUNDS} 轮",
                          f"- 性质窗口：MW±{WINDOW['mw']}, cLogP±{WINDOW['clogp']}, HBD±{WINDOW['hbd']}, "
                          f"HBA±{WINDOW['hba']}, RotB±{WINDOW['rotb']}, Rings±{WINDOW['rings']}",
                          f"- 抽样池可用分子数：{len(pool)}",
                          ""]

    for target, g in actives.groupby("target"):
        acts = g.reset_index(drop=True)
        act_fps = acts["fp"].tolist()
        n_act = len(acts)
        used = np.zeros(len(pool), dtype=bool)
        counts = np.zeros(n_act, dtype=int)
        assigned = []          # (active_idx, pool_idx)
        rng = np.random.default_rng(SEED)
        rounds_used = 0

        for rnd in range(MAX_ROUNDS):
            rounds_used = rnd + 1
            added = 0
            for ai in range(n_act):
                need = DECOYS_PER_ACTIVE - counts[ai]
                if need <= 0:
                    continue
                a = acts.iloc[ai]
                mask = ((~used)
                        & (np.abs(p_mw - a["mw"]) <= WINDOW["mw"])
                        & (np.abs(p_logp - a["clogp"]) <= WINDOW["clogp"])
                        & (np.abs(p_hbd - a["hbd"]) <= WINDOW["hbd"])
                        & (np.abs(p_hba - a["hba"]) <= WINDOW["hba"])
                        & (np.abs(p_rotb - a["rotb"]) <= WINDOW["rotb"])
                        & (np.abs(p_rings - a["rings"]) <= WINDOW["rings"]))
                cand = np.where(mask)[0]
                if len(cand) == 0:
                    continue
                # 与该靶点全部正样本的 Tanimoto 都要 < 0.4
                safe = []
                for ci in cand:
                    sims = DataStructs.BulkTanimotoSimilarity(pool_fps[ci], act_fps)
                    if max(sims) < TANIMOTO_MAX:
                        safe.append(ci)
                if not safe:
                    continue
                pick_n = min(need, len(safe))
                chosen = rng.choice(np.array(safe), size=pick_n, replace=False)
                for ci in chosen:
                    used[ci] = True
                    counts[ai] += 1
                    assigned.append((ai, int(ci)))
                    out_rows.append({"target": target,
                                     "molecule_chembl_id": pool.iloc[ci]["molecule_chembl_id"],
                                     "smiles": pool.iloc[ci]["smiles"]})
                    added += 1
            if counts.min() >= DECOYS_PER_ACTIVE:
                break
            if added == 0:
                break

        n_dec = len(assigned)
        ratio = n_dec / n_act if n_act else 0
        print(f"{target}: 正样本 {n_act}，诱饵 {n_dec}（{ratio:.1f}:1），轮次 {rounds_used}")
        logs.append(f"- {target}: 正样本 {n_act}，诱饵 {n_dec}（实际 {ratio:.1f}:1），"
                    f"抽样轮次 {rounds_used}"
                    + ("（未达 10:1，按实际记录）" if n_dec < n_act * DECOYS_PER_ACTIVE else "")
                    + ("；exploratory=True（正样本数不足，骨架切分后测试集过小）"
                       if target in EXPLORATORY else ""))

    dec = pd.DataFrame(out_rows)
    dec.to_csv(INTERIM / "decoys_raw.csv", index=False, encoding="utf-8-sig")

    logs += ["", "## 说明", "- 诱饵不跨靶点共享，各靶点独立构建；跨靶点允许同一分子重复，"
                 "pool_decoys.csv 每行只属于一个 target",
             "- 抽不满 10:1 的靶点按实际数量记录，未硬凑", ""]
    (LOGS / "decoy_log.md").write_text("\n".join(logs), encoding="utf-8")
    print(f"[10c] 诱饵总数 {len(dec)}")


if __name__ == "__main__":
    main()