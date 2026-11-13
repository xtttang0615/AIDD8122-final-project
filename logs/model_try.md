# 建模试错记录

- P4 初版 train CV AUC=1.0000 判定为特征选择泄漏 -> P4b 改嵌套 CV（固定池）。
- P4c 进一步发现固定池仍存在 pre-filtering bias（特征池来自全训练集 DEG）-> 改为完全嵌套 CV（折内重跑 limma），作为主结果。
- 分类器对比：logreg_l2=1.0000 vs random_forest=0.9957，沿用预设规则，未反复调整。
- 跨队列标准化只试了'各数据集内部逐基因 z-score'一种方式，未试其他标准化；结果照实报告。
- 未对 internal/external 做阈值或参数搜索；未改 C；未改 deg_* 文件。
- P3 火山图 adjustText 触发 matplotlib FancyArrowPatch 警告（自动改用 annotate），不修复仅记录。
