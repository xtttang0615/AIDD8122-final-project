# AIDD8122 Final Project

> This public repository contains the complete analysis code and result tables. The submitted PDF report is submitted separately; it carries the student's name and ID, which are omitted from this repository.

基于 GEO 外周血转录组数据，鉴定急性心肌梗死（AMI）诊断基因签名、按规则筛选可成药靶点，并用机器学习对候选分子打分排序。完整流水线见下文“脚本执行顺序”。

## 运行环境

### 解释器全路径（本机实测）

| 工具 | 路径（本机实测版本） |
|---|---|
| Python | `<Anaconda>\python.exe`（3.11） |
| Rscript | `<R install>\bin\Rscript.exe`（4.6.1） |
| git | `git`（PATH 中即可） |

> 提醒：`python` 命令会被 Windows 商店别名（App execution alias）拦截，脚本统一用解释器的完整路径调用。
> `src/07_signature.py` 通过环境变量 `RSCRIPT` 定位 Rscript（缺省回退到 PATH，再回退到上表路径），换机时设置 `RSCRIPT` 即可。
> Python 依赖见 `requirements.txt`：`pip install -r requirements.txt`；pandoc 与 R 环境需单独安装。

### R 包安装方式（供重装 / 换机参考）

```r
# 系统库 C:/R/R-4.6.1/library 不可写时，先装到用户库
.libPaths("C:/Users/<user>/AppData/Local/R/win-library/4.6")

# 清华镜像直连（不要用 BiocManager，其版本校验访问官方 bioconductor.org 会 SSL 失败）
options(repos = c(
  BioCsoft = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.23/bioc",
  BioCann  = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor/packages/3.23/data/annotation",
  CRAN     = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
))

# 默认 wininet 下载 https 会报 .inet_error，必须改用 libcurl
options(download.file.method = "libcurl")

install.packages(c("limma", "AnnotationDbi", "hgu133plus2.db"))
```

说明：

- Bioconductor 3.23 对应 R 4.6.1；安装时把 URL 里的 `3.23` 换成目标 R 版本对应的 Bioconductor 版本即可。
- `AnnotationDbi`、`hgu133plus2.db` 属于 Bioconductor 仓库，只设 CRAN 装不到，必须带上 `BioCsoft` / `BioCann`。
- 系统库不可写时，包会装到 `.libPaths()` 指定的用户库；R 启动会自动识别该用户库。

## 脚本执行顺序

| 阶段 | 脚本 | 说明 |
|---|---|---|
| 环境 | `00_check_env.py` / `00_check_env.R` | 生成 `logs/env_report.md` |
| P1 | `01_download.py` → `02_parse_geo.py` | 下载与解析 GEO（缓存 `data/raw/`） |
| P2 | `03_preprocess.py` → `04_qc_plot.py` | 探针→基因、质控图 |
| P3 | `05_deg_limma.R` → `06_deg_plots.py` | limma 差异分析 |
| P4 | `07_signature.py` → `08_roc.py`（`05b_limma_fold.R` 为折内 limma 辅助） | 签名建模与三层验证 |
| P5 | `09_target_selection.py` | 靶点规则化筛选 |
| P6 | `10a` → `10b` → `10c` → `10d` | 分子池构建与标准化 |
| P7 | `11_features.py` → `11_splits.py` → `12_ligand_model.py` → `13_ablation.py` → `14` → `15` → `16` | 配体基排序与评估 |
| P7b–P7d | `17_diagnostics.py`、`18b_hard_benchmark_v2.py`、`19_benchmark_normalize.py` | 泄漏诊断、硬负样本、可比性修正（旧版 `18_hard_benchmark.py` 已归档到 `archive/src/`） |
| P8 | `20` → `21` → `22` → `23` → `24` | 已上市药重定位与富集检验 |
| 出图 | `27_workflow_figure.py`、`28_supplements_v2.py`、`29_final_figures.py` | 流程图、补充图、图 8/9 定稿 |
| 审查 | `25_audit.py`、`30_verify_tables.py`、`31_verify_appendix_refs.py`、`32_audit_final.py`、`33_verify_references.py` | 分子身份与训练重叠核查、表格/附录/文献核验 |

## 数据来源与缓存

- GEO：GSE66360（发现集，CEC）、GSE48060（外部验证集，全血）、GPL570 注释。
- ChEMBL REST：活性、适应证、分子结构；响应缓存于 `data/raw/chembl/`（含 `approved/` 子目录），
  重跑不重复请求。缓存文件不入版本库。
- 结构化中间产物：`data/interim/*.parquet`、`data/processed/*.csv|*.npz`（均不入库，可由脚本重建）。

## 结果入口

- 转录组：`results/tables/deg_*.csv`、`model_performance.csv`、`signature_genes.csv`
- 靶点：`target_candidates.csv`（含 strict / relaxed 两套）
- 配体：`ligand_model_perf.csv`（主表）、`ligand_model_byfold.csv`、`ablation_summary.csv`、
  `ligand_model_hard_benchmark.csv`（含 legacy 归档）
- 重定位：`approved_pool.csv`、`approved_scores.csv`、`mi_enrichment_final.csv`、
  `repurposing_shortlist.csv`
- 诊断：`audit_*.csv`、`leak_check_diagnostics.csv`、`label_permutation_diagnostics.csv`、
  `cross_target_negatives.csv`
- 图：`results/figures/`（300 dpi）
- 修订说明：`docs/revision_record.md`、`docs/report_corrections.md`

## 项目定位

本项目的定位是**转录组靶点发现与配体药物排序流程的探索性研究及可靠性评价**，
不声称已验证的 AMI 治疗靶点，也不将排序结果解释为活性或结合概率预测。
详见 `docs/revision_record.md` 与 `docs/report_corrections.md`。