# release 自包含性验证

验证方式：把 `release/` 整目录复制到独立临时目录（`%TEMP%\aidd_release_check`），在副本内执行，
不触碰项目内的最终结果。

## 实际执行并通过

| 检查 | 结果 |
|---|---|
| `src/*.py` 语法编译（py_compile，38 个文件） | **0 失败** |
| `src/29_final_figures.py`（无网络、无训练） | 通过；输出 n=1823、median 0.020、p99 0.21468、max 0.89、>0.5 者 2，与交付表 `library_score_stats.csv` 一致 |
| `reports/build.py`（相对路径构建） | 通过；生成 `report_full.html` 3.12 MB，图片全部内联，未依赖本机绝对路径 |

## 静态检查通过（未执行）

| 检查 | 结果 |
|---|---|
| `src/28_supplements_v2.py` 输入 | 仅依赖 `results/tables/` 内已交付的 5 个 CSV，路径由 `parents[1]` 推导，不越出 release |
| `src/24_p8_figures.py` | 只输出 `fig7b_shortlist_heatmap.png`；fig9 已移交 29 号脚本，无重复输出 |
| `src/16_ligand_figures.py` | 依赖 `ligand_oof_predictions.csv`、`ligand_model_perf.csv`、`ablation_summary.csv`，均在 release 内 |
| 训练分支（01–23）的调用顺序与输入输出 | 按 README 顺序静态核对一致；`data/processed/` 已含两份固定切分 |

## 本轮明确未执行

- `01`–`23` 全部分析与训练步骤（按要求不重训）。
- R 脚本 `05_deg_limma.R`、`05b_limma_fold.R`（需要 R 与 Bioconductor 包；本机有 R 4.6.1，但本轮未运行）。
- 从原始 GEO/ChEMBL 数据开始的端到端复现。

## 真实阻碍

1. `11_features.py` 与后续配体步骤需要 `data/processed/ligand_features_ecfp4.npz`；该文件未随 release 发布，
   须由 `11_features.py` 从 `results/tables/pool_*.csv` 重建（无需网络）。
2. 需要 R 环境的两个脚本未纳入无 R 环境下的验证。
3. 未做「导入成功即端到端复现」的等价声明：以上仅为语法、依赖、路径与部分非训练步骤的实测。
