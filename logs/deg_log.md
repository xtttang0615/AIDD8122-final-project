# P3 limma 差异分析日志

- limma 版本：3.68.4
- train 样本数：43（case 21 / control 22）
- 主阈值：|logFC|>1 & adj.P.Val<0.05 -> 上调 370 / 下调 28 / 合计 398
- 放宽阈值：|logFC|>0.5 & adj.P.Val<0.05 -> 816
- internal_val 样本数：56
- internal 显著：600
- train/internal 重叠：190，重叠率 0.3808，方向一致率 1
- 基因名抽查：MARCH1, SEPT9；疑似被改写数：0
- top10 上调：SULF1, ITLN1, MAFB, PDK4, FAM198B, CSTA, PTX3, VCAN, NR4A3, EFEMP1
- top10 下调：COL5A3, NOG, ZNF578, WDR76, BTNL9, ANKEF1, CAMTA1, PAQR8, PTCH1, RORC


## 图与离群检查（06）
- 火山图：fig2g_volcano.png
- 热图前50：fig2h_heatmap_top50.png
- PCA top2000：fig2f_pca_top2000.png
- PC1 |z|>3 离群样本：[('GSM1620913', -3.09)]
