# -*- coding: utf-8 -*-
"""P6 共用工具：ChEMBL 请求（本地缓存 + 重试 + 请求日志）与 RDKit 标准化/描述符。

被 10a/10b/10c/10d 四个脚本导入，本身不产生结果文件。
"""
from pathlib import Path
import datetime
import hashlib
import json
import time

import requests
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, DataStructs, rdMolDescriptors
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize

ROOT = Path(__file__).resolve().parents[1]
CACHE_DIR = ROOT / "data" / "raw" / "chembl"
LOGS = ROOT / "logs"
REQ_LOG = LOGS / "api_requests.log"
INTERIM = ROOT / "data" / "interim"

BASE = "https://www.ebi.ac.uk/chembl/api/data"
TIMEOUT = 60
MAX_RETRIES = 2  # 首次 + 2 次重试，仍失败即抛出

CACHE_DIR.mkdir(parents=True, exist_ok=True)
INTERIM.mkdir(parents=True, exist_ok=True)
LOGS.mkdir(parents=True, exist_ok=True)


def _write_log(line):
    with open(REQ_LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def chembl_get(path, params, cache_name=None):
    """GET ChEMBL 接口；命中本地缓存则不发请求。返回 (json, from_cache)。"""
    if cache_name is None:
        key = hashlib.md5(f"{path}|{json.dumps(params, sort_keys=True)}".encode()).hexdigest()[:16]
        cache_name = f"{key}.json"
    cache_file = CACHE_DIR / cache_name
    if cache_file.exists():
        return json.loads(cache_file.read_text(encoding="utf-8")), True

    last = None
    for attempt in range(MAX_RETRIES + 1):
        t0 = time.time()
        try:
            r = requests.get(BASE + path, params=params, timeout=TIMEOUT,
                             headers={"User-Agent": "aidd8122-coursework"})
            dt = time.time() - t0
            if r.status_code == 200:
                data = r.json()
                cache_file.write_text(json.dumps(data), encoding="utf-8")
                _write_log(f"{datetime.datetime.now().isoformat(timespec='seconds')}\t200\t{dt:.2f}s\t"
                           f"{path}\tparams={json.dumps(params, sort_keys=True)[:180]}\tcache={cache_name}")
                return data, False
            last = f"HTTP {r.status_code}: {r.text[:200]}"
        except Exception as e:
            last = f"{type(e).__name__}: {e}"
        _write_log(f"{datetime.datetime.now().isoformat(timespec='seconds')}\tFAIL\t-\t{path}\t"
                   f"attempt={attempt+1}\t{last}")
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"{path} 请求失败：{last}")


def log_request_note(note):
    _write_log(f"{datetime.datetime.now().isoformat(timespec='seconds')}\tNOTE\t-\t{note}")


# ---------- RDKit ----------

METAL_Z = (set(range(21, 31)) | set(range(39, 49)) | set(range(57, 81)) | set(range(89, 104))
           | {3, 4, 11, 12, 13, 19, 20, 37, 38, 49, 50, 55, 56, 81, 82, 83, 84, 87, 88})

_LFC = rdMolStandardize.LargestFragmentChooser()
_UNCHARGER = rdMolStandardize.Uncharger()
_TAUTOMER = rdMolStandardize.TautomerEnumerator()
_FPGEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def standardize_mol(mol):
    """MolStandardize 流程：最大片段 -> 去盐 -> 电荷中和 -> canonical tautomer。

    只做结构标准化，不做分子性质过滤；失败返回 None。
    """
    if mol is None:
        return None
    try:
        mol = _LFC.choose(mol)
        mol = _UNCHARGER.uncharge(mol)
        mol = _TAUTOMER.Canonicalize(mol)
    except Exception:
        return None
    if mol is None or mol.GetNumAtoms() == 0:
        return None
    return mol


def prepare_smiles(smiles):
    """SMILES -> 标准化 mol，并给出丢弃原因（None 表示通过）。

    丢弃类别：empty_smiles / parse_fail / inorganic_metal / heavy_atoms_lt_8 / mw_gt_800
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return None, "empty_smiles"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "parse_fail"
    mol = standardize_mol(mol)
    if mol is None:
        return None, "parse_fail"
    if any(a.GetAtomicNum() in METAL_Z for a in mol.GetAtoms()):
        return None, "inorganic_metal"
    if mol.GetNumHeavyAtoms() < 8:
        return None, "heavy_atoms_lt_8"
    if Descriptors.MolWt(mol) > 800:
        return None, "mw_gt_800"
    return mol, None


def descriptors(mol):
    return {
        "smiles": Chem.MolToSmiles(mol),
        "inchikey14": Chem.MolToInchiKey(mol)[:14],
        "mw": round(Descriptors.MolWt(mol), 2),
        "clogp": round(Crippen.MolLogP(mol), 3),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
    }


def morgan_fp(mol):
    return _FPGEN.GetFingerprint(mol)


def tanimoto(fp, fps):
    return DataStructs.BulkTanimotoSimilarity(fp, fps)
