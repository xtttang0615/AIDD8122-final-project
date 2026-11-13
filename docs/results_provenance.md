# 报告结果对应表（数字 → 来源文件 → 筛选条件 → 分母）

所有数字均取自下列文件的最新版本；聊天记录只作核对线索。

| 报告数字 | 来源文件 | 筛选条件 / 分母 |
|---|---|---|
| 43 train（21/22）| `data/processed/GSE66360_meta.csv` | role=train |
| 56 internal（28/28）| 同上 | role=internal_val |
| 52 external（31/21）| `data/processed/GSE48060_meta.csv` | 全队列 |
| 398 DEG（370↑/28↓）| `results/tables/deg_train.csv` | \|logFC\|>1 & adj.P<0.05，train 内 |
| 190 / 600 显著（internal）| `results/tables/deg_internal.csv` | 同阈值，internal 内 |
| 重叠率 38.1% | `results/tables/deg_replication.csv` | 190 / mean(398, 600) |
| nested AUC 0.961（0.890–1.000）| `results/tables/model_performance.csv` | 行 set=train_nested_full，n=43 |
| fixed-pool AUC 0.950 | 同上 | 行 set=train_nested_fixedpool |
| optimistic AUC 1.000 | 同上 | 行 set=train_cv_optimistic（泄漏参考） |
| external AUC 0.608（0.442–0.770）| 同上 | 行 set=external，n=52 |
| cutoff 0.6001 | 同上 | 早期发现集 OOF 锁定的统一阈值 |
| 签名 10 基因 | `results/tables/signature_genes.csv` | 兜底 top-10（sel_freq≥50% 仅 1 个） |
| 靶点 28 个候选 / strict 1 / relaxed 7 | `results/tables/target_candidates*.csv` | OT 阈值；见 `logs/target_selection_log.md` |
| 正样本 800/273/244/44 | `results/tables/pool_actives.csv` | InChIKey14 去重；PDE4B cap 800 |
| 诱饵 6238/2658/2437/440 | `results/tables/pool_decoys.csv` | 逐靶点独立，1:1 抽样 |
| 配体 AP 0.9933/0.9679/0.9887 | `results/tables/ligand_model_perf.csv` | 行 model=rf, features=ecfp4（主模型） |
| 逐折波动 | `results/tables/ligand_model_byfold.csv` | 每折 1 行（PDK4 仅 1 次 hold-out） |
| 消融 A1/A2 | `results/tables/ablation_summary.csv` | A3 已标 no-op，`pr_auc` 为空 |
| 困难基准 AP | `results/tables/ligand_model_hard_benchmark.csv` | 新表 benchmark=hard_unified_folds；旧表 `_legacy` |
| MI 适应证漏斗 191→161→123→108→89 | `logs/p8_log.md`、`logs/mi_drugs_log.md` | ChEMBL drug_indication；InChIKey14 去重 |
| 富集检验（修订后）| `results/tables/mi_enrichment_final.csv` | 逐靶点排除训练重叠；PDE4B 86/1844, PFKFB3 89/1888, EPAS1 88/1879, consensus 86/1823 |
| 短名单 1822 候选 | `results/tables/repurposing_shortlist.csv`、`supp_S1_rank_consensus_top25.csv` | 排除 MI 适应证与任一模型训练重叠 |
| 控制药物 | `reports/supp_S3_control_molecules.csv` | 身份人工策展，分数/排名由脚本计算 |
| 得分-性质相关 | `reports/supp_S2_score_property_correlation.csv` | Spearman，n=1822 |
| 身份错配证据 | `results/tables/audit_mi_identity.csv` | 123 有结构 → 108 成功；截断差 14 个 ID |
| 训练重叠证据 | `results/tables/audit_overlap.csv` | InChIKey14；approved_pool 与负例重叠 61/17/26/6 |
