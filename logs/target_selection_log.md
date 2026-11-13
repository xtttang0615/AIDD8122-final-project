# 09 靶点规则化筛选日志

- 疾病实体：myocardial infarction（MONDO_0005068，由 Open Targets search 解析，未硬编码）
- 候选基因数（跨折签名池）：28
- 旧口径（OT overall >= 0.1 且 STRING degree >= 1）筛后：1
- strict（OT overall >= 0.1）筛后：1
- relaxed（OT overall >= 0.05 或 target_known_drugs_total > 0）筛后：7
- PDE4B 自检：target_known_drugs_total = 28（要求 > 0）

## 规则变更记录
- 原方案：hub = STRING degree 前 25% + ChEMBL >= 10 硬筛。
- 实际执行：简化为 OT overall >= 0.1 且 STRING degree >= 1（degree 门槛形同虚设，28 个基因中 27 个 degree 远大于 1），ChEMBL 仅作排序列。
- 变更理由：候选仅 28 个基因，hub 过滤区分度不足。
- P5b 追加：strict = OT >= 0.1；relaxed = OT >= 0.05 或 target_known_drugs_total > 0；known drugs 为靶点级字段（OT v4 drugAndClinicalCandidates.count，旧版 knownDrugs 已下线），ot_known_drug_disease 为 MI 疾病上下文内的已知药证据，两者语义不同，均已保留。

## 每接口调用量
- Open Targets：57 次（含 1 次疾病解析）
- STRING：28 次
- ChEMBL：42 次
- 失败记录：14 条（详见 logs/api_errors.md）

## 数据口径说明
- ChEMBL：每次取前 1000 条活性记录；chembl_n_molecules_first1000 为该批记录中去重的化合物数，chembl_n_activities_total 为接口返回的活性记录总数
- OT datatypeScores 只返回非零证据类型，字段为 NA 表示该证据类型无数据（零分省略）
- STRING：interaction_partners limit=1000，degree 为伙伴数（本次最大值远小于上限，未截断）

## 筛选结果（按 sel_freq -> OT overall -> max pChEMBL 降序）
| gene | sel_freq | OT overall | ot_known_drug_disease | target_known_drugs_total | STRING deg | max pChEMBL | strict | relaxed | relaxed_reason |
|---|---|---|---|---|---|---|---|---|---|
| PDE4B | 5 | 0.2655 |  | 28 | 121 | 11.0 | True | True | OT>=0.05 & known_drugs>0 |
| CXCL16 | 4 | 0.0626 |  | 0 | 133 |  | False | True | OT>=0.05 |
| GADD45A | 4 | 0.0089 |  | 0 | 128 |  | False | False |  |
| SKIL | 4 | 0.0042 |  | 0 | 67 |  | False | False |  |
| SYTL3 | 4 |  |  | 0 | 17 |  | False | False |  |
| EPAS1 | 3 | 0.0623 |  | 1 | 427 | 8.7 | False | True | OT>=0.05 & known_drugs>0 |
| METRNL | 2 | 0.0687 |  | 0 | 7 |  | False | True | OT>=0.05 |
| ZFP36 | 2 | 0.0413 |  | 0 | 195 |  | False | False |  |
| PDK4 | 2 | 0.0293 |  | 1 | 218 | 8.3 | False | True | known_drugs>0 |
| ITLN1 | 2 | 0.0292 |  | 0 | 51 |  | False | False |  |
| PFKFB3 | 1 | 0.0848 |  | 0 | 138 | 8.52 | False | True | OT>=0.05 |
| NR4A2 | 1 | 0.0446 |  | 0 | 131 | 9.1 | False | False |  |
| GZMB | 1 | 0.0384 |  | 0 | 424 | 9.3 | False | False |  |
| CH25H | 1 | 0.034 |  | 0 | 67 |  | False | False |  |
| CYP1B1 | 1 | 0.0222 |  | 0 | 160 | 10.37 | False | False |  |
| YBX3 | 1 | 0.0196 |  | 0 | 154 | 4.81 | False | False |  |
| AREG | 1 | 0.0187 |  | 0 | 269 | 6.46 | False | False |  |
| HBEGF | 1 | 0.0154 |  | 1 | 248 |  | False | True | known_drugs>0 |
| AKR1C3 | 1 | 0.009 |  | 0 | 122 | 9.47 | False | False |  |
| CD55 | 1 | 0.0075 |  | 0 | 151 |  | False | False |  |
| IRS2 | 1 | 0.006 |  | 0 | 196 |  | False | False |  |
| TNFAIP3 | 1 | 0.0058 |  | 0 | 295 |  | False | False |  |
| RGS1 | 1 | 0.0041 |  | 0 | 135 |  | False | False |  |
| SAMSN1 | 1 | 0.0033 |  | 0 | 38 | 8.19 | False | False |  |
| BCL9 | 1 |  |  | 0 | 108 | 5.8 | False | False |  |
| CHRM3-AS2 | 1 |  |  | 0 |  |  | False | False |  |
| SLITRK5 | 1 |  |  | 0 | 78 |  | False | False |  |
| SOX13 | 1 |  |  | 0 | 25 |  | False | False |  |
