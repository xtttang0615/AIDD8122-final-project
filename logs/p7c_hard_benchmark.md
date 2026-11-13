# P7c 硬负样本基准

硬负样本 = 现有随机诱饵 : 其他两靶点正样本 = 1:1；跨靶点分子按其骨架分配折号，骨架与本靶点正样本重叠者剔除
RF(500) + ECFP4，参数与 P7 完全一致，未调参；seed=42

## PDE4B
- 正样本 800；硬负样本 517（剔除骨架/InChIKey 冲突 0），随机诱饵抽样 517
- hard：PR-AUC=0.9976，ROC-AUC=0.9979，EF1=2.29，p@10=1.000
- easy：PR-AUC=0.9933；降幅 = -0.0043
- 正例率 = 0.4362；hard PR-AUC / 正例率 = 2.29x

- 该基准仍偏易（PR-AUC > 0.95），按要求不再加码

## PFKFB3
- 正样本 273；硬负样本 1044（剔除骨架/InChIKey 冲突 0），随机诱饵抽样 1044
- hard：PR-AUC=0.9689，ROC-AUC=0.9768，EF1=8.65，p@10=1.000
- easy：PR-AUC=0.9679；降幅 = -0.0010
- 正例率 = 0.1156；hard PR-AUC / 正例率 = 8.38x

- 该基准仍偏易（PR-AUC > 0.95），按要求不再加码

## EPAS1
- 正样本 244；硬负样本 1073（剔除骨架/InChIKey 冲突 0），随机诱饵抽样 1073
- hard：PR-AUC=0.9908，ROC-AUC=0.9966，EF1=9.80，p@10=1.000
- easy：PR-AUC=0.9887；降幅 = -0.0021
- 正例率 = 0.1021；hard PR-AUC / 正例率 = 9.70x

- 该基准仍偏易（PR-AUC > 0.95），按要求不再加码

## P7d 可比性修正

- 新增列：prevalence = n_pos/(n_pos+n_neg)；lift = PR-AUC / prevalence（随机排序的 PR-AUC 期望值即正例率）
- PDE4B easy 下采样对齐：从 6238 个负样本中随机抽 1034 个（seed=42），用已有 OOF 预测重算 20 次，PR-AUC = 0.9980 ± 0.0001

| target | benchmark | prevalence | PR-AUC | lift | ROC-AUC |
|---|---|---|---|---|---|
| PDE4B | hard | 0.4362 | 0.9976 | 2.29 | 0.9979 |
| PDE4B | easy | 0.1137 | 0.9933 | 8.74 | 0.9981 |
| PFKFB3 | hard | 0.1156 | 0.9689 | 8.38 | 0.9768 |
| PFKFB3 | easy | 0.0931 | 0.9679 | 10.39 | 0.9786 |
| EPAS1 | hard | 0.1021 | 0.9908 | 9.70 | 0.9966 |
| EPAS1 | easy | 0.0910 | 0.9887 | 10.86 | 0.9938 |
| PDE4B | easy_downsampled | 0.4362 | 0.9980 | 2.29 | 0.9980 |

ROC-AUC 对类别比例不敏感（0.9981 vs 0.9979），可作为跨基准比较的指标；PR-AUC 因正例率变化不可跨基准直接比较，须用 lift 或下采样对齐后比较。

### 结论（统一表述）

PFKFB3/EPAS1 正例率接近，hard 基准 PR-AUC 基本持平（lift 10.4→8.4、10.9→9.7），提示存在靶点特异判别；PDE4B 因正例率由 11.4% 升至 43.6%，其 PR-AUC 不可直接比较，下采样对齐后为 0.9980（lift 8.7→2.3）。

三分层排序：本靶点配体 > 其他靶点配体 > 随机小分子；两层基准均接近饱和，仅支持相对排序用途。
## P7c 修订（18b）：合并后统一骨架五折

### PDE4B
- 合并集合：1834（正 800 / 负 1034），骨架数 1032
- 统一 StratifiedGroupKFold(5, seed=42) 后跨折检查：scaffold 跨折 0，InChIKey14 跨折 0，标签冲突结构 0
- 剔除跨靶点分子（骨架/结构冲突）0
- hard_unified：PR-AUC=0.9975，ROC-AUC=0.9978，lift=2.29

### PFKFB3
- 合并集合：2361（正 273 / 负 2088），骨架数 1321
- 统一 StratifiedGroupKFold(5, seed=42) 后跨折检查：scaffold 跨折 0，InChIKey14 跨折 0，标签冲突结构 0
- 剔除跨靶点分子（骨架/结构冲突）0
- hard_unified：PR-AUC=0.9966，ROC-AUC=0.9992，lift=8.62

### EPAS1
- 合并集合：2390（正 244 / 负 2146），骨架数 1226
- 统一 StratifiedGroupKFold(5, seed=42) 后跨折检查：scaffold 跨折 0，InChIKey14 跨折 0，标签冲突结构 0
- 剔除跨靶点分子（骨架/结构冲突）0
- hard_unified：PR-AUC=0.9867，ROC-AUC=0.9941，lift=9.67

- 原切分（各靶点独立折号）保留在 data/processed/scaffold_splits.csv；新版切分写入 data/processed/hard_benchmark_splits.csv
- 其他靶点配体仅为“假定困难负例”，不等于该靶点实测非活性分子
- prevalence 变化会影响 AP，prevalence-matched 下采样仅作辅助描述
