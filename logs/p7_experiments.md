# P7 实验记录

## 2026-11-05T17:37:31 预注册主模型规则
- 默认主模型：XGBoost
- 规则：若三主靶点（PDE4B/PFKFB3/EPAS1）mean PR-AUC 上 XGB 与 RF 之差 < 0.02，取 RF（更简单者胜）
- 该规则在看结果之前写下，选定后不再更换

## 2026-11-05T17:39:11 主模型选定
- 三靶点 mean PR-AUC：XGB=0.9748，RF=0.9833，差值=0.0085（阈值 0.02）
- 结果：主模型 = **rf**

## 折间 PR-AUC（mean±SD，仅稳定性参考，不作为判据）
- PDE4B（rf）：0.9937 ± 0.0066（n_folds=5）
- PFKFB3（rf）：0.9521 ± 0.0769（n_folds=5）
- EPAS1（rf）：0.8935 ± 0.1924（n_folds=5）
- PDK4（rf）：1.0000 ± 0.0000（n_folds=1）

## 运行记录（时间/配置/结果/是否主结果）
- 2026-11-05T17:39:11 | EPAS1 | rf | ecfp4 | ROC=0.9938 PR=0.9887 | 主结果
- 2026-11-05T17:39:11 | EPAS1 | lr | ecfp4 | ROC=0.9926 PR=0.9869 | 对照
- 2026-11-05T17:39:11 | EPAS1 | xgb | ecfp4 | ROC=0.9919 PR=0.9739 | 对照
- 2026-11-05T17:39:11 | EPAS1 | lr | desc | ROC=0.8580 PR=0.2587 | 对照
- 2026-11-05T17:39:11 | EPAS1 | random | - | ROC=0.4932 PR=0.0877 | 对照
- 2026-11-05T17:39:11 | PDE4B | rf | ecfp4 | ROC=0.9981 PR=0.9933 | 主结果
- 2026-11-05T17:39:11 | PDE4B | lr | ecfp4 | ROC=0.9974 PR=0.9919 | 对照
- 2026-11-05T17:39:11 | PDE4B | xgb | ecfp4 | ROC=0.9970 PR=0.9893 | 对照
- 2026-11-05T17:39:11 | PDE4B | lr | desc | ROC=0.6416 PR=0.2261 | 对照
- 2026-11-05T17:39:11 | PDE4B | random | - | ROC=0.4995 PR=0.1135 | 对照
- 2026-11-05T17:39:11 | PDK4 | rf | ecfp4 | ROC=1.0000 PR=1.0000 | 主结果
- 2026-11-05T17:39:11 | PDK4 | lr | ecfp4 | ROC=0.9988 PR=0.9889 | 对照
- 2026-11-05T17:39:11 | PDK4 | xgb | ecfp4 | ROC=0.9737 PR=0.9211 | 对照
- 2026-11-05T17:39:11 | PDK4 | lr | desc | ROC=0.6547 PR=0.1477 | 对照
- 2026-11-05T17:39:11 | PDK4 | random | - | ROC=0.5451 PR=0.1091 | 对照
- 2026-11-05T17:39:11 | PFKFB3 | rf | ecfp4 | ROC=0.9786 PR=0.9679 | 主结果
- 2026-11-05T17:39:11 | PFKFB3 | lr | ecfp4 | ROC=0.9759 PR=0.9624 | 对照
- 2026-11-05T17:39:11 | PFKFB3 | xgb | ecfp4 | ROC=0.9689 PR=0.9611 | 对照
- 2026-11-05T17:39:11 | PFKFB3 | lr | desc | ROC=0.7927 PR=0.2962 | 对照
- 2026-11-05T17:39:11 | PFKFB3 | random | - | ROC=0.4825 PR=0.0874 | 对照

## 2026-11-05T17:41:18 消融（13）
- 主模型：rf；切分复用 scaffold_splits.csv，未重新生成
- A1_features | PDE4B | rf | ecfp4 | default | PR=0.9933 ROC=0.9981 EF1=8.80
- A1_features | PFKFB3 | rf | ecfp4 | default | PR=0.9679 ROC=0.9786 EF1=10.74
- A1_features | EPAS1 | rf | ecfp4 | default | PR=0.9887 ROC=0.9938 EF1=10.99
- A1_features | PDE4B | rf | desc | default | PR=0.6757 ROC=0.8875 EF1=8.80
- A1_features | PFKFB3 | rf | desc | default | PR=0.5513 ROC=0.8697 EF1=10.37
- A1_features | EPAS1 | rf | desc | default | PR=0.6523 ROC=0.9374 EF1=8.95
- A1_features | PDE4B | rf | both | default | PR=0.9933 ROC=0.9982 EF1=8.80
- A1_features | PFKFB3 | rf | both | default | PR=0.9672 ROC=0.9785 EF1=10.74
- A1_features | EPAS1 | rf | both | default | PR=0.9893 ROC=0.9957 EF1=10.99
- A2_model | PDE4B | lr | ecfp4 | default | PR=0.9919 ROC=0.9974 EF1=8.80
- A2_model | PFKFB3 | lr | ecfp4 | default | PR=0.9624 ROC=0.9759 EF1=10.74
- A2_model | EPAS1 | lr | ecfp4 | default | PR=0.9869 ROC=0.9926 EF1=10.99
- A2_model | PDE4B | rf | ecfp4 | default | PR=0.9933 ROC=0.9981 EF1=8.80
- A2_model | PFKFB3 | rf | ecfp4 | default | PR=0.9679 ROC=0.9786 EF1=10.74
- A2_model | EPAS1 | rf | ecfp4 | default | PR=0.9887 ROC=0.9938 EF1=10.99
- A2_model | PDE4B | xgb | ecfp4 | default | PR=0.9893 ROC=0.9970 EF1=8.80
- A2_model | PFKFB3 | xgb | ecfp4 | default | PR=0.9611 ROC=0.9689 EF1=10.74
- A2_model | EPAS1 | xgb | ecfp4 | default | PR=0.9739 ROC=0.9919 EF1=10.99
- A3_threshold | PDE4B | rf | ecfp4 | pchembl>=7 | PR=nan ROC=nan EF1=nan
- A3 note：仅比较同一切分下的相对变化，正负比已改变，不作绝对性能解读

