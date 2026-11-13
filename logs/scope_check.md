# 口径核对

- 打分表 1994；MI 标记 89；非 MI 1905；与训练数据重叠 82；**共同候选 1823**。
- 与 `repurposing_shortlist.csv` 的候选差集：0 个。
- 合并性质后缺失 mw 的分子：0（缺者即为图注需说明的分子）。
- 排名定义：PDE4B 排名在共同候选集合内按 PDE4B 得分降序，并列取同分同名次（method="min"），分母 1823。
- 聚合：主要用 mean percentile（三靶点百分位均值），Borda（秩和）仅作补充比较。

                       comparison  spearman_rho  top20_overlap  top50_overlap
     mean_percentile vs prob_mean        0.9519             11             30
mean_percentile vs borda_rank_sum        0.9999             20             49
