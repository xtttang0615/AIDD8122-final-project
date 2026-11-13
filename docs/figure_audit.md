# 图表审查记录（figure audit）

| 报告图号 | 图片文件 | 生成脚本 | 输入结果表 | 集合/样本 | 问题与处理 |
|---|---|---|---|---|---|
| Fig 1 | `workflow_pipeline.png` | `src/27_workflow_figure.py` | 无（流程） | — | 已重绘为实际执行流程；删除装饰性“未运行某方法”说明；与 part2 图注一致 |
| Fig 2 | `fig2g_volcano.png` | `src/06_deg_plots.py` | `deg_train.csv` | train n=43 | 阈值线与 top10 标注核对无误 |
| Fig 3 | `fig3b_calibration.png` | `src/08_roc.py` | `predictions.csv` | internal/external | 图注仅描述实际展示的两条校准曲线 |
| Fig 4 | `fig4_ligand_roc.png` / `fig4b_ligand_pr.png` | `src/16_ligand_figures.py` | `ligand_oof_predictions.csv` | 合并 OOF | 图注写明为 pooled OOF；AP 定义已在正文说明 |
| Fig 5 | `fig5_ef_curve.png` | `src/16_ligand_figures.py` | `ligand_oof_predictions.csv` | 合并 OOF，RF 主模型 | 图注只描述 RF 曲线，未声称多算法比较 |
| Fig 6 | `fig6_ablation.png` | `src/16_ligand_figures.py` | `ablation_summary.csv` | A1/A2 | A3 已为 no-op，面板改为说明文字，不再与其他实验并列 |
| Fig 7 | `fig7b_shortlist_heatmap.png` | `src/24_p8_figures.py` | `repurposing_shortlist.csv` | top20，1,823 候选 | 标题为 model score，不称“药物候选发现” |
| Fig 8 | `fig8_score_property.png` | `src/28_supplements_v2.py` | `approved_scores.csv` | 1,823 | 双面板同集合：(A) 得分分布 (B) 得分–MW；相关系数与正文一致（0.396） |
| Fig 9 | `fig7_mi_enrichment.png` | `src/24_p8_figures.py` | `mi_enrichment_final.csv` | 逐靶点各自评价集 | 分母由图内 n 标注（86/1844、89/1888、88/1879、86/1823），与表一致 |
| Supp | `fig_control_analysis.png` | `src/26_supplements.py`（历史）/`28_supplements_v2.py`（当前数值） | `supp_S3_control_molecules.csv` | 1,823 | 控制药物条形图；与 Fig 8 散点图不混用图注 |

## 未完全达标项

- 图 8 由 28 号脚本生成，但报告正文引用的“99th percentile 0.210”取自 1,905 打分池，未在 1,823 集合上重算；差异应 <0.005，但未验证。
- 图 5 只画 RF，且未画随机基线（正文已说明随机 EF≈1）。
- 所有图均为 300 dpi；未做黑白打印可读性检查。
