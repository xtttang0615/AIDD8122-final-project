# P4 建模日志（P4c：完全嵌套 CV + 标准化敏感性）

- 随机种子：SEED=42（外层折 rng=SEED+fold_index）
- 全训练集特征池：398 个 DEG（deg_train.csv 中 sig==TRUE）
- 主筛选 bootstrap：200 次；嵌套各折 bootstrap：100 次
- 筛选规则：RF(n_estimators=500) top-20 ∩ LASSO(L1,liblinear,C=1.0) 非零集，sel_freq >= 50%；不足 5 个按 sel_freq -> rf_importance_mean 降序取 top-10
- sel_freq >= 50% 的基因总数：1
- 最终签名（全训练集）：PDE4B, GADD45A, CXCL16, ZFP36, SKIL, SYTL3, SAMSN1, EPAS1, ITLN1, RGS1
- 分类器：logreg_l2=1.0000 vs random_forest=0.9957 -> 选 logreg_l2（未调 C）

## 完全嵌套 CV（主结果）
- train_nested_full AUC = 0.9610 [0.8896, 1.0000]
- 每折流程：折内训练部分 limma(|logFC|>1 & adj.P<0.05) -> bootstrap 100 次稳定性筛选 -> L2 LR

- fold 1：折内 DEG 236 个，签名 10 个（兜底 top-10）；基因：PDE4B, SAMSN1, GADD45A, ITLN1, SLITRK5, CHRM3-AS2, CXCL16, EPAS1, PDK4, SKIL
- fold 2：折内 DEG 296 个，签名 10 个（兜底 top-10）；基因：PDE4B, CXCL16, GADD45A, AREG, SYTL3, METRNL, SKIL, CD55, EPAS1, HBEGF
- fold 3：折内 DEG 263 个，签名 10 个（兜底 top-10）；基因：PDE4B, SYTL3, RGS1, GADD45A, NR4A2, CH25H, ITLN1, BCL9, CXCL16, CYP1B1
- fold 4：折内 DEG 442 个，签名 10 个（兜底 top-10）；基因：CXCL16, PDE4B, SKIL, METRNL, ZFP36, SYTL3, GZMB, AKR1C3, PFKFB3, PDK4
- fold 5：折内 DEG 270 个，签名 10 个（兜底 top-10）；基因：PDE4B, GADD45A, EPAS1, SYTL3, ZFP36, SOX13, YBX3, IRS2, TNFAIP3, SKIL

## 对照结果
- train_nested_fixedpool AUC = 0.9502 [0.8723, 1.0000]：特征池由全训练集 DEG 定义，存在 pre-filtering bias，不作主结果
- train_cv_optimistic AUC = 1.0000：全量筛选 + 同批样本 CV，存在选择泄漏，不作主结果

## 三层验证与阈值
- 统一阈值 cutoff=0.6001（P4 锁定，三层共用），另并列 cutoff=0.5；未对 internal/external 单独寻阈值
- external 缺失基因数：0（缺失用 train 中位数填补）

- train_nested_full: n=43 (case 21/control 22), AUC=0.961 [0.8896, 1.0], @0.6001 sens=0.9048 spec=0.9545 acc=0.9302, @0.5 sens=0.9524 spec=0.9091 acc=0.9302 | 主结果（完全嵌套 CV）
- train_nested_fixedpool: n=43 (case 21/control 22), AUC=0.9502 [0.8723, 1.0], @0.6001 sens=0.9048 spec=0.9091 acc=0.907, @0.5 sens=0.9524 spec=0.9091 acc=0.9302 | 特征池由全训练集 DEG 定义，存在 pre-filtering bias
- train_cv_optimistic: n=43 (case 21/control 22), AUC=1.0 [1.0, 1.0], @0.6001 sens=0.9524 spec=1.0 acc=0.9767, @0.5 sens=1.0 spec=1.0 acc=1.0 | 存在选择泄漏，不作为主结果
- internal_val: n=56 (case 28/control 28), AUC=0.8418 [0.7194, 0.949], @0.6001 sens=0.6429 spec=1.0 acc=0.8214, @0.5 sens=0.6429 spec=0.9643 acc=0.8036
- external: n=52 (case 31/control 21), AUC=0.6083 [0.4424, 0.7696], @0.6001 sens=0.0 spec=1.0 acc=0.4038, @0.5 sens=0.0 spec=1.0 acc=0.4038
- external_zscored: n=52 (case 31/control 21), AUC=0.5745 [0.4132, 0.742], @0.6001 sens=0.5806 spec=0.619 acc=0.5962, @0.5 sens=0.6129 spec=0.619 acc=0.6154 | post-hoc 敏感性分析，逐数据集内部 z-score，未使用任何标签

## 跨队列标准化敏感性分析（post-hoc）
- 做法：GSE66360 用 train 43 例的逐基因均值/标准差、GSE48060 用其自身 52 例的逐基因均值/标准差，各自标准化后训练与预测；未使用任何标签信息。
- external AUC：0.6083 -> 0.5745；mean_prob：0.0079 -> 0.5178

## 校准
- internal_val: slope=0.2987, intercept=0.9087, Brier=0.1606, mean_prob=0.3454
- external: slope=0.4048, intercept=2.6108, Brier=0.5849, mean_prob=0.0079
- external_zscored: slope=0.0889, intercept=0.3686, Brier=0.3145, mean_prob=0.5178

## 结论（照实记录）
- bootstrap 200 次中 sel_freq >= 50% 仅 1 个基因（PDE4B: 0.785），签名稳定性不足。
- external AUC = 0.6083（95% CI 0.4424-0.7696），跨细胞类型（CEC -> 全血）迁移性能有限。

## 图（08）
- fig3_roc.png：四条 ROC 同图（train_nested_full / internal_val / external / external_zscored），
  标 AUC 与 95% CI，对角虚线为 chance
- fig3b_calibration.png：校准曲线，external 为主，另附 internal_val 与 Brier 分数

## 附1：兜底规则的精确定义与触发原因

**规则定义（pick_signature）**
1. 先取 sel_freq >= 0.50 的基因集合 S（sel_freq = 200 次 bootstrap 中
   "RF top-20 ∩ LASSO 非零集"命中的比例）。
2. 若 |S| >= 5，则 S 全部入选，签名内部按 sel_freq 降序排列。
3. 若 |S| < 5（本次 |S| = 1），触发兜底：不做频率门槛过滤，直接把特征池全部基因按
   第一排序键 sel_freq 降序、第二排序键 rf_importance_mean 降序（`np.lexsort((-rf_imp, -sel_freq))`）
   排序，取前 10 个作为签名。等价于"取稳定性排名最靠前的 10 个基因"。

**触发原因**：200 次 bootstrap 中只有 PDE4B 的 sel_freq 达到 0.785，其余基因全部 < 0.5
（第 2 名 GADD45A 仅 0.430）。说明在 43 例样本上，RF top-20 与 LASSO 非零集的交集高度不稳定，
单一 50% 频率门槛无法产生 >= 5 个基因的签名，故按预设规则走 top-10 兜底。

**五折完全嵌套 CV 中同样每折都触发兜底**：折内 DEG 236–442 个，但 sel_freq >= 50% 的基因均不足 5 个。

## 附2：五折签名 10 基因两两 overlap 矩阵（上三角，交集数/Jaccard）

五折签名合计涉及 28 个不同基因，任意两折交集 3–5 个，Jaccard 0.176–0.333。

| | F1 | F2 | F3 | F4 | F5 |
|---|---|---|---|---|---|
| F1 | 1.00 | 5 / 15 = 0.333 | 4 / 16 = 0.250 | 4 / 16 = 0.250 | 4 / 16 = 0.250 |
| F2 | | 1.00 | 4 / 16 = 0.250 | 5 / 15 = 0.333 | 5 / 15 = 0.333 |
| F3 | | | 1.00 | 3 / 17 = 0.176 | 3 / 17 = 0.176 |
| F4 | | | | 1.00 | 4 / 16 = 0.250 |
| F5 | | | | | 1.00 |

交集数矩阵（对称，仅列上三角的原始计数）：F1-F2=5，F1-F3=4，F1-F4=4，F1-F5=4，F2-F3=4，
F2-F4=5，F2-F5=5，F3-F4=3，F3-F5=3，F4-F5=4。

跨折出现次数：PDE4B 5/5；GADD45A、CXCL16、SKIL、SYTL3 各 4/5；EPAS1 3/5；
ITLN1、PDK4、METRNL、ZFP36 各 2/5；其余 17 个基因各 1/5。
即五折之间仅 PDE4B 完全一致，签名集合层面的重合度中等偏低（Jaccard <= 0.333）。

## 附3：external_zscored 的统计解读（照实记录）

1. **95% CI (0.4132 - 0.7420) 包含 0.5**，即跨数据集（CEC -> 全血）判别与随机不可区分；
   点估计 0.5745 不能作为有效迁移的证据。
2. **Brier = 0.3145 > 0.25**（0.25 相当于"所有样本恒定预测 0.5"的平凡基线），
   说明标准化后的概率预测精度不如平凡基线；**校准斜率 0.0889 约等于 0**，
   预测概率几乎不携带可分信息。
3. 因此：external 与 external_zscored 的概率输出**不可用于临床阈值决策**，
   仅 AUC 可作为排序能力的参考指标。原 external（未标准化）mean_prob = 0.0079、
   Brier = 0.5849，问题更严重（概率整体左移导致 cutoff 完全失效）。

