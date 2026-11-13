# -*- coding: utf-8 -*-
"""11 特征工程：ECFP4(2048) + 9 维描述符，一次性缓存。

输入：results/tables/pool_actives.csv、pool_decoys.csv
输出：data/processed/ligand_features_ecfp4.npz、data/processed/ligand_features_desc.csv
运行：C:/Anaconda/Anaconda3/python.exe src/11_features.py

特征矩阵只建一次，12/13/14/15/16 全部复用，禁止各脚本重算。
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
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors, rdFingerprintGenerator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ligand_utils import DESC_CSV, NPZ, PROC, ROOT, TAB  # noqa: E402
from pool_utils import prepare_smiles  # noqa: E402

RDLogger.DisableLog("rdApp.*")

FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def main():
    act = pd.read_csv(TAB / "pool_actives.csv", encoding="utf-8-sig")
    dec = pd.read_csv(TAB / "pool_decoys.csv", encoding="utf-8-sig")
    rows = pd.concat([act[["molecule_chembl_id", "smiles"]],
                      dec[["molecule_chembl_id", "smiles"]]], ignore_index=True)
    rows = rows.drop_duplicates("molecule_chembl_id", keep="first")
    print(f"待计算分子：{len(rows)}（正样本+诱饵去重后）")

    ids, fps, descs = [], [], []
    bad = 0
    for r in rows.itertuples(index=False):
        mol, _ = prepare_smiles(r.smiles)
        if mol is None:
            bad += 1
            continue
        ids.append(r.molecule_chembl_id)
        fps.append(FPGEN.GetFingerprintAsNumPy(mol).astype(np.uint8))
        descs.append({
            "molecule_chembl_id": r.molecule_chembl_id,
            "mw": Descriptors.MolWt(mol),
            "clogp": Crippen.MolLogP(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "tpsa": rdMolDescriptors.CalcTPSA(mol),
            "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "rings": rdMolDescriptors.CalcNumRings(mol),
            "arom_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
            "fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        })

    X = np.vstack(fps)
    PROC.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(NPZ, ids=np.array(ids), ecfp4=X)
    pd.DataFrame(descs).to_csv(DESC_CSV, index=False, encoding="utf-8-sig")
    print(f"特征矩阵：{X.shape[0]} x {X.shape[1]}，无法标准化 {bad}")
    print("[11_features] 完成")


if __name__ == "__main__":
    main()