# 09 靶点规则化筛选规格（09_target_selection.py）

## 1. 候选范围

- 输入基因集合 = **跨折签名池 28 个基因**（完全嵌套 CV 五折签名 + 最终签名去重后的并集），
  不是只取最终 10 个基因。理由：签名跨折不稳定，仅用固定折的 10 基因做靶点筛选是人为截断。
- 每个基因附两列标注：
  - `sel_freq`：该基因在 5 个外层折签名中出现的次数（1–5）
  - `in_final_signature`：是否属于全训练集最终签名 10 基因（TRUE/FALSE）
- **筛选不按频率硬删**，频率信息全部保留进排序，供人工判断。

## 2. 数据源（三个都用，全部 REST，免费无密钥）

| 数据源 | 端点 | 取用字段 |
|---|---|---|
| Open Targets Platform | `https://api.platform.opentargets.org/api/v4/graphql` | overall association score；`datatypeScores` 分解（重点 `genetic_association`、`known_drug`）；`tractability`（可成药性注释）；`proteinIds`（UniProt，用于 ChEMBL 映射） |
| STRING | `https://string-db.org/api/json/interaction_partners` | `degree`（互作伙伴数）、`combined_score`（最高值） |
| ChEMBL | `https://www.ebi.ac.uk/chembl/api/data` | 活性化合物数（去重 `molecule_chembl_id`）、`max pChEMBL` |

- 疾病实体（心肌梗死）通过 Open Targets search 接口解析得到 EFO/MONDO id，**不硬编码、不凭记忆写**。
- 每请求 60s 超时、3 次指数退避重试。
- 映射失败或接口失败的基因：写入 `logs/api_errors.md`（基因名 + 端点 + 错误原文），
  对应字段留 `NA`，**绝不编造**，也不使用记忆中的靶点-药物关系（红线 2/3 继续生效）。

## 3. 硬筛选规则（取宽）

- 入选条件：`OT overall >= 0.1` **且** `STRING degree >= 1`。
- ChEMBL 活性化合物数**不做硬筛**（好靶点常无化学起点），仅作排序列。
- 可成药家族（kinase / GPCR / ion channel / protease 等）用 Open Targets `tractability`
  字段作**注释列**，不做硬筛。
- 理由：候选池仅 28 个基因，硬筛过狠会清空列表。

## 4. 排序与标注

- `target_candidates.csv` 排序键（均降序）：
  1. `sel_freq`（跨折出现次数）
  2. `ot_overall_score`
  3. `chembl_max_pchembl`
- `known_drug` 标注取自 Open Targets `datatypeScores` 中 id 为 `known_drug` 的字段，属数据来源，
  不使用任何模型记忆。

## 5. 产出

- `results/tables/target_candidates.csv`：全部列（含 `NA` 标注），含被硬筛剔除的记录，
  另加 `pass_filter` 列标明是否入选。
- `logs/target_selection_log.md`：筛前/筛后基因数、每接口调用量、失败清单、阈值说明。
- 运行结束后打印 Top 10 全文供人工审。

## 6. 约束

- 不改 `data/processed/**`、不改 `results/tables/deg_*`、不改特征池与阈值。
- 不安装新依赖（仅用 requests + pandas）。
