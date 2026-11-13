# P6 分子池构建日志

## 一、正样本漏斗（每靶点）

| target | exploratory | 原始分子 | 空 SMILES | 标准化后 | InChIKey14 去重 | 最终正样本 | capped |
|---|---|---|---|---|---|---|---|
| EPAS1 | False | 273 | 0 | 270 | 244 | 244 | 否 |
| PDE4B | False | 1720 | 8 | 1703 | 1651 | 800 | 是（capped at 800） |
| PDK4 | True | 44 | 0 | 44 | 44 | 44 | 否 |
| PFKFB3 | False | 286 | 0 | 286 | 273 | 273 | 否 |

- cap 明细（PDE4B）：1720 → InChIKey14 去重 1651 → 取 pChEMBL 前 800 → 800

- exploratory 靶点：PDK4（理由：正样本数不足，骨架切分后测试集过小）

## 二、丢弃原因计数

- 正样本：{'mw_gt_800': 12, 'empty_smiles': 8}
- 诱饵：无
- MI 金标准：{'mw_gt_800': 7, 'heavy_atoms_lt_8': 7, 'empty_smiles': 38, 'inorganic_metal': 1}
- empty_smiles（接口未返回结构）与 parse_fail（RDKit 无法解析）分开计数

## 三、诱饵

| target | 正样本 | 诱饵数 | 比例 |
|---|---|---|---|
| PDE4B | 800 | 6238 | 7.8:1 |
| PFKFB3 | 273 | 2658 | 9.74:1 |
| EPAS1 | 244 | 2437 | 9.99:1 |
| PDK4 | 44 | 440 | 10.0:1 |

## 四、MI 药物金标准

- 漏斗：161 个药物 → 有结构 123 → 标准化后可用 108
- 无结构 38（结构缺失者排除出评估集，保留在名单内备查）
- 查询方式：efo_term__contains=myocardial infarction
- max_phase：approved(4)=138，clinical(1-3)=23

### efo_term 分布（按药物数）

- myocardial infarction: 135
- acute myocardial infarction: 55
- anterolateral myocardial infarction: 1

## 五、泄漏自检（评估集 vs 训练正样本，InChIKey14）

- EPAS1: 正样本 244，与 MI 金标准重叠 0
- PDE4B: 正样本 800，与 MI 金标准重叠 0
- PDK4: 正样本 44，与 MI 金标准重叠 0
- PFKFB3: 正样本 273，与 MI 金标准重叠 0
- **评估集与训练集重叠数 = 0，P7 将剔除**
