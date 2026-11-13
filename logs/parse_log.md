# P1 解析日志

## GSE66360

### 数据集概况

- 样本数：99
- 样本来源（!Sample_source_name_ch1 唯一值）：['CD146+ Circulating Endothelial Cells']
- 样本标题（!Sample_title 唯一值）：['Muse_1_DISCOVERY_Control', 'Muse_2_DISCOVERY_Control', 'Muse_3_DISCOVERY_Control', 'Muse_4_DISCOVERY_Control', 'Muse_5_DISCOVERY_Control', 'Muse_6_DISCOVERY_Control', 'Muse_7_DISCOVERY_Control', 'Muse_8_DISCOVERY_Control', 'Muse_9_DISCOVERY_Control', 'Muse_10_DISCOVERY_Control', 'Muse_11_DISCOVERY_Control', 'Muse_12_DISCOVERY_Control', 'Muse_13_DISCOVERY_Control', 'Muse_14_DISCOVERY_Myocardial Infarction', 'Muse_15_DISCOVERY_Myocardial Infarction', 'Muse_16_DISCOVERY_Myocardial Infarction', 'Muse_17_DISCOVERY_Myocardial Infarction', 'Muse_18_DISCOVERY_Myocardial Infarction', 'Muse_19_DISCOVERY_Myocardial Infarction', 'Muse_20_DISCOVERY_Myocardial Infarction', 'Muse_21_DISCOVERY_Myocardial Infarction', 'Muse_22_DISCOVERY_Myocardial Infarction', 'Muse_23_DISCOVERY_Myocardial Infarction', 'Muse_24_DISCOVERY_Myocardial Infarction', 'Muse_25_DISCOVERY_Myocardial Infarction', 'Muse_26_DISCOVERY_Myocardial Infarction', 'Muse_27_DISCOVERY_Myocardial Infarction', 'Muse_28_DISCOVERY_Myocardial Infarction', 'Muse_29_DISCOVERY_Myocardial Infarction', 'Muse_30_DISCOVERY_Myocardial Infarction', 'Muse_31_DISCOVERY_Myocardial Infarction', 'Muse_32_DISCOVERY_Myocardial Infarction', 'Muse_33_DISCOVERY_Myocardial Infarction', 'Muse_34_DISCOVERY_Myocardial Infarction', 'Muse_35_DISCOVERY_Control', 'Muse_36_DISCOVERY_Control', 'Muse_37_DISCOVERY_Control', 'Muse_38_DISCOVERY_Control', 'Muse_39_DISCOVERY_Control', 'Muse_40_DISCOVERY_Control', 'Muse_41_DISCOVERY_Control', 'Muse_42_DISCOVERY_Control', 'Muse_43_DISCOVERY_Control', 'Muse_44_VALIDATION_Control', 'Muse_45_VALIDATION_Control', 'Muse_46_VALIDATION_Control', 'Muse_47_VALIDATION_Control', 'Muse_48_VALIDATION_Control', 'Muse_49_VALIDATION_Control', 'Muse_50_VALIDATION_Control', 'Muse_51_VALIDATION_Control', 'Muse_52_VALIDATION_Control', 'Muse_53_VALIDATION_Control', 'Muse_54_VALIDATION_Control', 'Muse_55_VALIDATION_Control', 'Muse_56_VALIDATION_Control', 'Muse_57_VALIDATION_Control', 'Muse_58_VALIDATION_Control', 'Muse_59_VALIDATION_Myocardial Infarction', 'Muse_60_VALIDATION_Myocardial Infarction', 'Muse_61_VALIDATION_Myocardial Infarction', 'Muse_62_VALIDATION_Myocardial Infarction', 'Muse_63_VALIDATION_Myocardial Infarction', 'Muse_64_VALIDATION_Myocardial Infarction', 'Muse_65_VALIDATION_Myocardial Infarction', 'Muse_66_VALIDATION_Myocardial Infarction', 'Muse_67_VALIDATION_Myocardial Infarction', 'Muse_68_VALIDATION_Myocardial Infarction', 'Muse_69_VALIDATION_Myocardial Infarction', 'Muse_70_VALIDATION_Myocardial Infarction', 'Muse_71_VALIDATION_Myocardial Infarction', 'Muse_72_VALIDATION_Myocardial Infarction', 'Muse_73_VALIDATION_Myocardial Infarction', 'Muse_74_VALIDATION_Control', 'Muse_75_VALIDATION_Control', 'Muse_76_VALIDATION_Control', 'Muse_77_VALIDATION_Control', 'Muse_78_VALIDATION_Control', 'Muse_79_VALIDATION_Control', 'Muse_80_VALIDATION_Control', 'Muse_81_VALIDATION_Control', 'Muse_82_VALIDATION_Control', 'Muse_83_VALIDATION_Control', 'Muse_84_VALIDATION_Control', 'Muse_85_VALIDATION_Control', 'Muse_86_VALIDATION_Control', 'Muse_87_VALIDATION_Myocardial Infarction', 'Muse_88_VALIDATION_Myocardial Infarction', 'Muse_89_VALIDATION_Myocardial Infarction', 'Muse_90_VALIDATION_Myocardial Infarction', 'Muse_91_VALIDATION_Myocardial Infarction', 'Muse_92_VALIDATION_Myocardial Infarction', 'Muse_93_VALIDATION_Myocardial Infarction', 'Muse_94_VALIDATION_Myocardial Infarction', 'Muse_95_VALIDATION_Myocardial Infarction', 'Muse_96_VALIDATION_Myocardial Infarction', 'Muse_97_VALIDATION_Myocardial Infarction', 'Muse_98_VALIDATION_Myocardial Infarction', 'Muse_99_VALIDATION_Myocardial Infarction']

### 分组映射

- 规则：case <- `acute myocardial infarction` / `stemi` / `ami`；control <- `healthy` / `normal` / `control`
- control: 50
- case: 49
- 无法判断样本：无

### 尺度探测

- min=1.4904  median=5.9729  max=14.9000  1%=2.3886  99%=12.3561
- 判定：已 log2 转换（14.9000 <=100）
- 说明：本阶段不做任何变换，仅记录。

### characteristics 键值频次

- 键 `cell type`（99 次）
  - `CD146+ Circulating Endothelial Cells` -> 99

### discovery/validation 检测

- characteristics 中未找到 discovery/validation 字段。
- **!Sample_title 中含官方 DISCOVERY/VALIDATION 标注**：{'VALIDATION': np.int64(56), 'DISCOVERY': np.int64(43)}
- 按既定约定，split 列仍按数据集划分：GSE66360=discovery，GSE48060=validation。

## GSE48060

### 数据集概况

- 样本数：52
- 样本来源（!Sample_source_name_ch1 唯一值）：['Peripheral blood, patient with recurrent events', 'Peripheral blood, patient without recurrent events', 'Peripheral blood, normal control']
- 样本标题（!Sample_title 唯一值）：['Peripheral blood, patient with recurrent events [N21]', 'Peripheral blood, patient without recurrent events [N23]', 'Peripheral blood, patient without recurrent events [N24]', 'Peripheral blood, patient without recurrent events [N25]', 'Peripheral blood, patient without recurrent events [N26]', 'Peripheral blood, patient without recurrent events [N27]', 'Peripheral blood, patient with recurrent events [N28]', 'Peripheral blood, patient without recurrent events [N29]', 'Peripheral blood, patient without recurrent events [N30]', 'Peripheral blood, patient with recurrent events [N31]', 'Peripheral blood, patient without recurrent events [N32]', 'Peripheral blood, patient without recurrent events [N33]', 'Peripheral blood, patient without recurrent events [N34]', 'Peripheral blood, patient without recurrent events [N35]', 'Peripheral blood, patient without recurrent events [N36]', 'Peripheral blood, patient without recurrent events [N37]', 'Peripheral blood, patient with recurrent events [N38]', 'Peripheral blood, patient without recurrent events [N39]', 'Peripheral blood, patient without recurrent events [N40]', 'Peripheral blood, patient without recurrent events [N41]', 'Peripheral blood, patient without recurrent events [N42]', 'Peripheral blood, patient without recurrent events [N43]', 'Peripheral blood, patient without recurrent events [N44]', 'Peripheral blood, patient without recurrent events [N45]', 'Peripheral blood, patient without recurrent events [N46]', 'Peripheral blood, patient with recurrent events [N47]', 'Peripheral blood, patient without recurrent events [N48]', 'Peripheral blood, patient without recurrent events [N49]', 'Peripheral blood, patient without recurrent events [N50]', 'Peripheral blood, patient without recurrent events [N51]', 'Peripheral blood, normal control [N52]', 'Peripheral blood, normal control [N53]', 'Peripheral blood, normal control [N54]', 'Peripheral blood, normal control [N55]', 'Peripheral blood, normal control [N56]', 'Peripheral blood, normal control [N57]', 'Peripheral blood, normal control [N58]', 'Peripheral blood, normal control [N59]', 'Peripheral blood, normal control [N61]', 'Peripheral blood, normal control [N62]', 'Peripheral blood, normal control [N63]', 'Peripheral blood, normal control [N64]', 'Peripheral blood, normal control [N65]', 'Peripheral blood, normal control [N66]', 'Peripheral blood, normal control [N67]', 'Peripheral blood, normal control [N68]', 'Peripheral blood, normal control [N69]', 'Peripheral blood, normal control [N70]', 'Peripheral blood, normal control [N71]', 'Peripheral blood, normal control [N73]', 'Peripheral blood, patient without recurrent events [N74]', 'Peripheral blood, normal control [N75]']

### 分组映射

- 规则：case <- `acute myocardial infarction` / `stemi` / `ami`；control <- `healthy` / `normal` / `control`
- case: 31
- control: 21
- 无法判断样本：无

### 尺度探测

- min=2.1341  median=5.4088  max=15.4071  1%=2.6290  99%=11.4995
- 判定：已 log2 转换（15.4071 <=100）
- 说明：本阶段不做任何变换，仅记录。

### characteristics 键值频次

- 键 `disease status`（52 次）
  - `patient without recurrent events` -> 26
  - `normal control` -> 21
  - `patient with recurrent events` -> 5

### discovery/validation 检测

- characteristics 中未找到 discovery/validation 字段。
- split 按数据集划分：GSE66360=discovery，GSE48060=validation（外部验证集，不参与调参）。

## 探针合并规则（预留）

- 本阶段（P1）不执行多探针合并。
- 后续 03_preprocess.py 将实现：多探针映射同一基因时取 **median**。
- 03 脚本已预留 `merge_probes()` 实现位置（P2 提供）。

## 异常字段说明

- 本阶段未做过滤，矩阵中存在探针缺失/非数值时以 NaN 记录，待 P2 处理。
