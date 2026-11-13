# 10a 活性数据拉取日志

- PDE4B (CHEMBL275): 过滤后活性记录 2202（接口 total_count=2202），去重分子 1720，其中有 canonical_smiles 1712
- CXCL16: no_ligand_data（无 ChEMBL target id）
- EPAS1 (CHEMBL1744522): 过滤后活性记录 363（接口 total_count=363），去重分子 273，其中有 canonical_smiles 273
- METRNL: no_ligand_data（无 ChEMBL target id）
- PDK4 (CHEMBL4141): 过滤后活性记录 52（接口 total_count=52），去重分子 44，其中有 canonical_smiles 44
- PFKFB3 (CHEMBL2331053): 过滤后活性记录 348（接口 total_count=348），去重分子 286，其中有 canonical_smiles 286
- HBEGF (CHEMBL3286070): no_ligand_data（0 条满足 pChEMBL>=6 的活性）

## 无配体靶点（no_ligand_data，不进建模）
- CXCL16, METRNL, HBEGF（有 target id CHEMBL3286070，但 0 条合格活性）
