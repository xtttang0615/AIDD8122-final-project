# P7b 诊断日志

## 一、集合与切分体检

                        check target  pre_standardization  post_standardization detail  n_cross_fold_keys
    actives_vs_decoys_overlap  EPAS1                  0.0                   0.0                       NaN
    actives_vs_decoys_overlap  PDE4B                  0.0                   0.0                       NaN
    actives_vs_decoys_overlap   PDK4                  0.0                   0.0                       NaN
    actives_vs_decoys_overlap PFKFB3                  0.0                   0.0                       NaN
cross_fold_scaffold/group_key  PDE4B                  NaN                   NaN                       0.0
        cross_fold_inchikey14  PDE4B                  NaN                   NaN                       0.0
cross_fold_scaffold/group_key PFKFB3                  NaN                   NaN                       0.0
        cross_fold_inchikey14 PFKFB3                  NaN                   NaN                       0.0
cross_fold_scaffold/group_key  EPAS1                  NaN                   NaN                       0.0
        cross_fold_inchikey14  EPAS1                  NaN                   NaN                       0.0
cross_fold_scaffold/group_key   PDK4                  NaN                   NaN                       0.0
        cross_fold_inchikey14   PDK4                  NaN                   NaN                       0.0

## 二、标签置换对照

target  roc_auc  pr_auc  prevalence    n
 PDE4B   0.4988  0.1303      0.1137 7038
PFKFB3   0.4914  0.0944      0.0931 2931
 EPAS1   0.3859  0.0687      0.0910 2681

## 三、消融全表

    ablation target model features    setting  pr_auc  roc_auc     ef1                                                   note
 A1_features  PDE4B    rf    ecfp4    default  0.9933   0.9981  8.7975                                                    NaN
 A1_features PFKFB3    rf    ecfp4    default  0.9679   0.9786 10.7363                                                    NaN
 A1_features  EPAS1    rf    ecfp4    default  0.9887   0.9938 10.9877                                                    NaN
 A1_features  PDE4B    rf     desc    default  0.6757   0.8875  8.7975                                                    NaN
 A1_features PFKFB3    rf     desc    default  0.5513   0.8697 10.3660                                                    NaN
 A1_features  EPAS1    rf     desc    default  0.6523   0.9374  8.9529                                                    NaN
 A1_features  PDE4B    rf     both    default  0.9933   0.9982  8.7975                                                    NaN
 A1_features PFKFB3    rf     both    default  0.9672   0.9785 10.7363                                                    NaN
 A1_features  EPAS1    rf     both    default  0.9893   0.9957 10.9877                                                    NaN
    A2_model  PDE4B    lr    ecfp4    default  0.9919   0.9974  8.7975                                                    NaN
    A2_model PFKFB3    lr    ecfp4    default  0.9624   0.9759 10.7363                                                    NaN
    A2_model  EPAS1    lr    ecfp4    default  0.9869   0.9926 10.9877                                                    NaN
    A2_model  PDE4B    rf    ecfp4    default  0.9933   0.9981  8.7975                                                    NaN
    A2_model PFKFB3    rf    ecfp4    default  0.9679   0.9786 10.7363                                                    NaN
    A2_model  EPAS1    rf    ecfp4    default  0.9887   0.9938 10.9877                                                    NaN
    A2_model  PDE4B   xgb    ecfp4    default  0.9893   0.9970  8.7975                                                    NaN
    A2_model PFKFB3   xgb    ecfp4    default  0.9611   0.9689 10.7363                                                    NaN
    A2_model  EPAS1   xgb    ecfp4    default  0.9739   0.9919 10.9877                                                    NaN
A3_threshold  PDE4B    rf    ecfp4 pchembl>=7  0.9933   0.9981  8.7975 仅比较同一切分下的相对变化，正负比 800:6238 -> 800:6238，正负比已改变，不作绝对性能解读

### 单维描述符 AUC（auc_flipped = 方向归一后的可分性）

target descriptor    auc  auc_flipped
PFKFB3 arom_rings 0.7155       0.7155
 EPAS1       rotb 0.3219       0.6781
 EPAS1      rings 0.3247       0.6753
 EPAS1 arom_rings 0.3348       0.6652
 PDE4B        hba 0.6354       0.6354

## 四、跨靶点硬负样本

model_target donor_target  n_negatives  frac_pred_pos  median_percentile  cutoff
       PDE4B       PFKFB3          273         0.1758            88.5337   0.210
       PDE4B        EPAS1          244         0.1598            88.6189   0.210
      PFKFB3        PDE4B          800         0.0412            90.0375   0.114

## 五、Top 候选新颖性

- novel 个数：{'PDE4B': 20, 'PFKFB3': 20, 'EPAS1': 20}

## 六、随机基线（100 次置换 mean±SD）

target    n  prevalence  pr_mean  pr_sd  ef1_mean  ef1_sd
 PDE4B 7038      0.1137   0.1154 0.0033    1.0205  0.3254
PFKFB3 2931      0.0931   0.0955 0.0065    0.9515  0.6751
 EPAS1 2681      0.0910   0.0941 0.0063    1.0377  0.5679
