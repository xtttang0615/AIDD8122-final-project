# P8 已上市药物重定位

## 一、候选池漏斗

- ChEMBL max_phase=4 小分子总数：3475
- 拉取记录：3475
- 无结构：164
- 标准化失败（无法解析/金属/重原子<8/MW>800）：228
- 标准化后：3083
- InChIKey14 去重后：1988
- 与训练正样本（pool_actives）重叠剔除：3
- 与 MI 金标准重叠剔除：80
- **最终候选池：1905**
- 打分覆盖：1994/1994 个分子三靶点均可打分（100.0%）
- 打分表含候选池 1905 个 + MI 金标准 89 个（后者仅供富集检验与短名单排除）
- 打分模型：RF(500)+ECFP4，与 P7 主模型参数一致，未调参

## 三、MI 药物富集检验（定稿）

- 池内 MI 药 89 个，其余已上市药 1905 个
   target  n_MI  n_nonMI  median_MI  median_nonMI    U_p  auc_like    ef1    ef5 conclusion
    PDE4B    89     1905     0.0180        0.0200 0.5943    0.4925 2.2404 0.6721     未观察到富集
   PFKFB3    89     1905     0.0100        0.0120 0.7027    0.4834 1.1202 0.4481     未观察到富集
    EPAS1    89     1905     0.0100        0.0100 0.4104    0.5071 1.1202 0.6721     未观察到富集
consensus    89     1905     0.0173        0.0173 0.7442    0.4794 0.0000 1.1202     未观察到富集

## 四、重定位短名单（top 20 预览，已排除 MI 金标准）

| rank | pref_name | chembl_id | prob_mean |
|---|---|---|---|
| 1 | CRISABOROLE | CHEMBL484785 | 0.3907 |
| 2 | TAVABOROLE | CHEMBL443052 | 0.2060 |
| 3 | LENVATINIB | CHEMBL1289601 | 0.1947 |
| 4 | PEFICITINIB | CHEMBL3137308 | 0.1560 |
| 5 | TUCATINIB | CHEMBL3989868 | 0.1527 |
| 6 | REVUMENIB | CHEMBL4650827 | 0.1380 |
| 7 | SORAFENIB | CHEMBL1336 | 0.1373 |
| 8 | PONATINIB | CHEMBL1171837 | 0.1353 |
| 9 | DORAVIRINE | CHEMBL2364608 | 0.1340 |
| 10 | TIVOZANIB | CHEMBL1289494 | 0.1313 |
| 11 | CAPMATINIB | CHEMBL3188267 | 0.1300 |
| 12 | IVOSIDENIB | CHEMBL3989958 | 0.1287 |
| 13 | REGORAFENIB | CHEMBL1946170 | 0.1267 |
| 14 | CABOZANTINIB S-MALATE | CHEMBL2103868 | 0.1240 |
| 15 | SIMEPREVIR | CHEMBL501849 | 0.1227 |
| 16 | DASATINIB ANHYDROUS | CHEMBL1421 | 0.1207 |
| 17 | TELAPREVIR | CHEMBL231813 | 0.1153 |
| 18 | ETRAVIRINE | CHEMBL308954 | 0.1140 |
| 19 | ATOGEPANT | CHEMBL3991065 | 0.1133 |
| 20 | FINERENONE | CHEMBL2181927 | 0.1127 |

机制/适应症/文献列一律留空，待人工核对。

> 模型饱和，排序仅用于相对优先级，不作活性预测。

## 五、化学型驱动的假阳性（TAVABOROLE，更正）

- 此前将 TAVABOROLE 描述为“硼酸类 PDE4 抑制剂”，该表述**错误**，现更正如下。
- **TAVABOROLE 的确认靶点为亮氨酰-tRNA 合成酶（LeuRS，抗真菌药），非 PDE4；
  其 PDE4B 高分源于与 crisaborole 共享硼酸化学型，属化学型驱动的假阳性，
  恰好说明排序反映结构相似性而非确认药理学。**

## 六、阳性对照外部验证

- **roflumilast（CHEMBL193240）**：InChIKey14 = `MNDBXUUTURYVHR`，pChEMBL = 9.82，
  命中 PDE4B 训练集（pool_actives）。
- **apremilast（CHEMBL514800）**：InChIKey14 = `IMOZEMNVLZVGJZ`，pChEMBL = 7.78，
  同样命中 PDE4B 训练集（pool_actives）。
- 两者按预注册规则在候选池构建阶段即被剔除（与训练正样本重叠），
  **目的是防止循环论证**——若让已知训练集分子参与排序，排序首位必然是它自己。
- **crisaborole（CHEMBL484785）** 不在训练集（其 InChIKey14 未出现在 pool_actives 中），
  模型对其独立打分：`prob_PDE4B = 0.89`，在候选池 1905 个分子中排名 **1/1905**
  （prob_mean = 0.3907，prob_min_rank = 1）。
- 这构成“未见过该分子的模型召回了已上市 PDE4 抑制剂”，是**本研究唯一的阳性对照证据**；
  由于只有 1 个对照分子，只能作为可行性提示，不足以支撑泛化性能结论。

## 七、MI 药计数口径统一

- 全文出现的“108”指标准化后可用数；按 InChIKey14 去重后为 **89**
  （同一药物的不同盐型/重复记录被合并）。
- `mi_enrichment_final.csv` 统一使用 89；引用时应表述为“108 个可用 → 去重后 89 个”，
  两者不矛盾。
## 三（修订）、MI 适应证关联检查（排除训练数据重叠后）

- 集合定义：ChEMBL 中具有 MI 适应证记录的分子（非“已获批 AMI 药物金标准”）
- PDE4B: 剔除与该模型训练数据重叠 64 个分子，剩余评价集 1930（MI 86 / 非 MI 1844）
- PFKFB3: 剔除与该模型训练数据重叠 17 个分子，剩余评价集 1977（MI 89 / 非 MI 1888）
- EPAS1: 剔除与该模型训练数据重叠 27 个分子，剩余评价集 1967（MI 88 / 非 MI 1879）
- consensus: 剔除与该模型训练数据重叠 85 个分子，剩余评价集 1909（MI 86 / 非 MI 1823）

   target  n_MI  n_nonMI  median_MI  median_nonMI    U_p  auc_like    ef1    ef5 conclusion
    PDE4B    86     1844     0.0180         0.020 0.5886    0.4929 2.3623 0.4675     未观察到富集
   PFKFB3    89     1888     0.0100         0.012 0.7358    0.4803 1.1107 0.6731     未观察到富集
    EPAS1    88     1879     0.0100         0.010 0.4830    0.5013 1.1176 0.6843     未观察到富集
consensus    86     1823     0.0173         0.018 0.7495    0.4786 0.0000 1.1683     未观察到富集

## 四、重定位短名单（top 20 预览，已排除 MI 金标准）

| rank | pref_name | chembl_id | prob_mean |
|---|---|---|---|
| 1 | CRISABOROLE | CHEMBL484785 | 0.3907 |
| 2 | TAVABOROLE | CHEMBL443052 | 0.2060 |
| 3 | LENVATINIB | CHEMBL1289601 | 0.1947 |
| 4 | PEFICITINIB | CHEMBL3137308 | 0.1560 |
| 5 | TUCATINIB | CHEMBL3989868 | 0.1527 |
| 6 | REVUMENIB | CHEMBL4650827 | 0.1380 |
| 7 | SORAFENIB | CHEMBL1336 | 0.1373 |
| 8 | PONATINIB | CHEMBL1171837 | 0.1353 |
| 9 | DORAVIRINE | CHEMBL2364608 | 0.1340 |
| 10 | TIVOZANIB | CHEMBL1289494 | 0.1313 |
| 11 | CAPMATINIB | CHEMBL3188267 | 0.1300 |
| 12 | IVOSIDENIB | CHEMBL3989958 | 0.1287 |
| 13 | REGORAFENIB | CHEMBL1946170 | 0.1267 |
| 14 | CABOZANTINIB S-MALATE | CHEMBL2103868 | 0.1240 |
| 15 | SIMEPREVIR | CHEMBL501849 | 0.1227 |
| 16 | DASATINIB ANHYDROUS | CHEMBL1421 | 0.1207 |
| 17 | TELAPREVIR | CHEMBL231813 | 0.1153 |
| 18 | ETRAVIRINE | CHEMBL308954 | 0.1140 |
| 19 | ATOGEPANT | CHEMBL3991065 | 0.1133 |
| 20 | FINERENONE | CHEMBL2181927 | 0.1127 |

机制/适应症/文献列一律留空，待人工核对。

> 模型饱和，排序仅用于相对优先级，不作活性预测。

## 十、审查后口径修订（P9）

- 集合命名：本集合统一称“ChEMBL 中具有 MI 适应证记录的分子”，不再称“已获批 AMI 药物金标准”。
  检索方式 efo_term__contains="myocardial infarction"（efo_id=EFO_0000612 命中 0 条）。
  分子总体 max_phase 不等同 AMI 适应证获批。
- 完整漏斗：191 条适应证记录 -> 161 分子 -> 有结构 123 -> 标准化可用 108 -> 去重 89。
- 评价独立性：原候选池仅排除训练正例（pool_actives），未排除训练负例（pool_decoys）。
  核查显示 approved_pool 与训练负例重叠 PDE4B 61 / PFKFB3 17 / EPAS1 26 / PDK4 6。
  现逐靶点剔除与该模型全部训练数据重叠的分子，共识分析使用三模型均未见过的共同集合。
- 富集检验定位为探索性适应证关联检查；四个维度均未观察到富集，结论不变。
- 短名单聚合排名同步排除训练重叠分子（1905 -> 1823），前 20 排序未变。
