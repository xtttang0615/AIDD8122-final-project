# 未来提交文件清单（release manifest）

原则：只含最终代码、最终结果、报告生成源与必要说明；旧脚本与 legacy 输出本地归档，不进提交包。

## 1. 代码（最终运行链）

| 顺序 | 文件 | 用途 |
|---|---|---|
| 0 | `src/00_check_env.py` / `.R` | 环境自检 |
| 1 | `src/01_download.py` → `src/02_parse_geo.py` | GEO 下载与解析 |
| 2 | `src/03_preprocess.py` → `src/04_qc_plot.py` | 探针映射、质控图 |
| 3 | `src/05_deg_limma.R` → `src/06_deg_plots.py` | limma DE 与图 |
| 4 | `src/07_signature.py`（+ `src/05b_limma_fold.R`）→ `src/08_roc.py` | 签名与三层验证 |
| 5 | `src/09_target_selection.py` | 靶点筛选 |
| 6 | `src/10a` → `10b` → `10c` → `10d` | 分子池构建与标准化 |
| 7 | `src/11_features.py` → `src/11_splits.py` → `src/12_ligand_model.py` → `src/13_ablation.py` → `14` → `15` → `16` | 配体基准 |
| 8 | `src/17_diagnostics.py` | 泄漏/偏置诊断 |
| 9 | `src/18b_hard_benchmark_v2.py`（**唯一入口**）→ `src/19_benchmark_normalize.py` | 硬负样本与可比性 |
| 10 | `src/20` → `21` → `22` → `23` → `24` | 重定位与富集 |
| 11 | `src/25_audit.py` → `src/26_supplements.py` → `src/28_supplements_v2.py` | 审查与补充表（28 为当前口径） |
| 12 | `src/27_workflow_figure.py` | 流程图 |
| — | `src/pool_utils.py`、`src/ligand_utils.py` | 共用工具 |

## 2. 归档（不进提交包，保留在本机）

- `src/18_hard_benchmark.py`（被 18b 取代：合并后未统一骨架分折）
- `src/26_supplements.py` 初版数值（被 28 取代：Borda 方向错误、候选集含 PDK4 排除）
- `results/tables/mi_drug_scores_legacy_misaligned.csv`（身份错配版）
- `results/tables/ligand_model_hard_benchmark_legacy.csv`（旧折号）
- `reports/AIDD8122_Final_Report.docx/.pdf`、`AIDD8122_Final_Report_revised.docx/.pdf`、`reports/_build/`
- 根目录 `student ID+name+final project.docx`（保留为原件，不覆盖）

## 3. 数据

- `data/raw/`：GEO 原始文件与 ChEMBL 缓存（不入库）；由 01/10a/10b/20 脚本按 URL 重建。
- `data/processed/scaffold_splits.csv`、`hard_benchmark_splits.csv`：**需入库**（固定切分，复现依赖）。
- `data/processed/ligand_features_ecfp4.npz`：可由 `11_features.py` 重建，不入库。

## 4. 结果与报告

- 结果表：`results/tables/*.csv`（剔除 `*_legacy*`、`*_legacy_misaligned*`）
- 图：`results/figures/*.png`、`workflow_pipeline.png`、`reports/fig8_score_property.png`
- 报告源：`reports/part1.html`–`part4.html`、`reports/build.py`
- 最终交付：`reports/AIDD8122_Final_Submission.docx/.pdf`
- 说明文档：`README.md`、`docs/*.md`

## 5. 明确排除

`__pycache__/`、`*.pyc`、`.screenshots/`、`final project/`（个人材料子目录）、
`15天每日行动计划.docx`、`TRAE指令手册.md`、课程 PDF/PPTX、`reports/_cand_desc.pkl`、
`reports/_build/`。
