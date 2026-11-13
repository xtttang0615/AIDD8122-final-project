# 修订前核查记录（P9-audit）

## A. 14 号脚本分子身份对应

                                                                                                                                                                                                                                      0
n_with_smiles                                                                                                                                                                                                                       123
n_feature_ok                                                                                                                                                                                                                        108
n_failed                                                                                                                                                                                                                             15
identity_mismatch_by_truncation                                                                                                                                                                                                    True
n_ids_only_in_truncation                                                                                                                                                                                                             14
n_ids_only_in_success                                                                                                                                                                                                                14
failed_ids                         CHEMBL1174;CHEMBL117822;CHEMBL1200574;CHEMBL1200644;CHEMBL1200689;CHEMBL1200731;CHEMBL1234886;CHEMBL160;CHEMBL1644695;CHEMBL395429;CHEMBL3989857;CHEMBL412531;CHEMBL414357;CHEMBL4802129;CHEMBL93268
failed_names                     EPTIFIBATIDE;ETHYLENE;SODIUM CHLORIDE;FONDAPARINUX SODIUM;NITRIC OXIDE;POTASSIUM CHLORIDE;OXYGEN;CYCLOSPORINE;SODIUM IODIDE;OXYTOCIN;MAGNESIUM SULFATE;THYMOSIN;EXENATIDE;BALIXAFORTIDE;SODIUM NITRITE

失败分子：
molecule_chembl_id           pref_name           reason
        CHEMBL1174        EPTIFIBATIDE        mw_gt_800
      CHEMBL117822            ETHYLENE heavy_atoms_lt_8
     CHEMBL1200574     SODIUM CHLORIDE heavy_atoms_lt_8
     CHEMBL1200644 FONDAPARINUX SODIUM        mw_gt_800
     CHEMBL1200689        NITRIC OXIDE heavy_atoms_lt_8
     CHEMBL1200731  POTASSIUM CHLORIDE  inorganic_metal
     CHEMBL1234886              OXYGEN heavy_atoms_lt_8
         CHEMBL160        CYCLOSPORINE        mw_gt_800
     CHEMBL1644695       SODIUM IODIDE heavy_atoms_lt_8
      CHEMBL395429            OXYTOCIN        mw_gt_800
     CHEMBL3989857   MAGNESIUM SULFATE heavy_atoms_lt_8
      CHEMBL412531            THYMOSIN        mw_gt_800
      CHEMBL414357           EXENATIDE        mw_gt_800
     CHEMBL4802129       BALIXAFORTIDE        mw_gt_800
       CHEMBL93268      SODIUM NITRITE heavy_atoms_lt_8

## C. 候选/评价集与训练数据重叠（InChIKey14）

            set target  n_set  n_overlap_train_pos  n_overlap_train_neg  n_overlap_all_train                                                                                                                                                                                                                                                                                                 overlap_ids
  approved_pool  PDE4B   1905                    0                   61                   61 ADEBPBSSDYVVLD;AGHANLSBXUWXTB;BCFGMOOMADDAQU;BJWJSFCBOIUATF;BRUQQQPBMZOVGD;BWWVAEOLVKTZFQ;CARRCRVKWDEWCD;CBQGYUDMJHNJBX;CHNUOJQWGUIOLD;CKHSUZOTJQWRKZ;DGBIGWXXNGSACT;FABPRXSRWADJSP;GPMXUUPHFNMNDH;GZUITABIAKMVPG;HRWCVUIFMSZDJS;HSUGRBWQSSZJOP;IEDVJHCEMCRBQM;IMZMKUWMOSJXDT;JFPVXVDWJQMJEE;KJONHKAYOJNZEC
  approved_pool PFKFB3   1905                    0                   17                   17                                              ADEBPBSSDYVVLD;BJWJSFCBOIUATF;CGIGDMFJXJATDK;GPMXUUPHFNMNDH;GZUITABIAKMVPG;IENZQIKPVFGBNW;JZCPYUJPEARBJL;MJIHNNLFOKEZEW;MUTNCGKQJGXKEM;NBQKINXMPLXUET;OVCDSSHSILBFBN;PCZOHLXUXFIOCF;UHTHHESEBZOYNR;WDLWHQDACQUCJR;XCPGHVQEEXUHNC;XMAYWYJOQHXEEK;XMXHEBAFVSFQEX
  approved_pool  EPAS1   1905                    0                   26                   26 BCEHBSKCWLPMDN;CARRCRVKWDEWCD;CGIGDMFJXJATDK;CHNUOJQWGUIOLD;DGBIGWXXNGSACT;GDLIGKIOYRNHDA;GPMXUUPHFNMNDH;HSUGRBWQSSZJOP;JOFWLTCLBGQGBO;LKJPYSCBVHEWIU;MBUVEWMHONZEQD;MJIHNNLFOKEZEW;MNJVRJDLRVPLFE;OUBCNLGXQFSTLU;PCZOHLXUXFIOCF;RAKDWPCGCPYIRC;RGCVKNLCSQQDEP;RKLNONIVDFXQRX;RZEKVGVHFLEQIL;URKOMYMAXPYINW
  approved_pool   PDK4   1905                    0                    6                    6                                                                                                                                                                                                                   LKJPYSCBVHEWIU;MJIHNNLFOKEZEW;MNJVRJDLRVPLFE;QAGYKUNXZHXKMR;QYIYFLOTGYLRGG;YFAGHNZHGGCZAX
approved_scores  PDE4B   1994                    0                   64                   64 ADEBPBSSDYVVLD;AGHANLSBXUWXTB;AHOUBRCZNHFOSL;BCFGMOOMADDAQU;BJWJSFCBOIUATF;BRUQQQPBMZOVGD;BWWVAEOLVKTZFQ;CARRCRVKWDEWCD;CBQGYUDMJHNJBX;CHNUOJQWGUIOLD;CKHSUZOTJQWRKZ;DGBIGWXXNGSACT;FABPRXSRWADJSP;GPMXUUPHFNMNDH;GZUITABIAKMVPG;HRWCVUIFMSZDJS;HSUGRBWQSSZJOP;IAKHMKGGTNLKSZ;IEDVJHCEMCRBQM;IMZMKUWMOSJXDT
approved_scores PFKFB3   1994                    0                   17                   17                                              ADEBPBSSDYVVLD;BJWJSFCBOIUATF;CGIGDMFJXJATDK;GPMXUUPHFNMNDH;GZUITABIAKMVPG;IENZQIKPVFGBNW;JZCPYUJPEARBJL;MJIHNNLFOKEZEW;MUTNCGKQJGXKEM;NBQKINXMPLXUET;OVCDSSHSILBFBN;PCZOHLXUXFIOCF;UHTHHESEBZOYNR;WDLWHQDACQUCJR;XCPGHVQEEXUHNC;XMAYWYJOQHXEEK;XMXHEBAFVSFQEX
approved_scores  EPAS1   1994                    0                   27                   27 AHOUBRCZNHFOSL;BCEHBSKCWLPMDN;CARRCRVKWDEWCD;CGIGDMFJXJATDK;CHNUOJQWGUIOLD;DGBIGWXXNGSACT;GDLIGKIOYRNHDA;GPMXUUPHFNMNDH;HSUGRBWQSSZJOP;JOFWLTCLBGQGBO;LKJPYSCBVHEWIU;MBUVEWMHONZEQD;MJIHNNLFOKEZEW;MNJVRJDLRVPLFE;OUBCNLGXQFSTLU;PCZOHLXUXFIOCF;RAKDWPCGCPYIRC;RGCVKNLCSQQDEP;RKLNONIVDFXQRX;RZEKVGVHFLEQIL
approved_scores   PDK4   1994                    0                    6                    6                                                                                                                                                                                                                   LKJPYSCBVHEWIU;MJIHNNLFOKEZEW;MNJVRJDLRVPLFE;QAGYKUNXZHXKMR;QYIYFLOTGYLRGG;YFAGHNZHGGCZAX
        mi_gold  PDE4B     89                    0                    3                    3                                                                                                                                                                                                                                                                AHOUBRCZNHFOSL;IAKHMKGGTNLKSZ;RUZYUOTYCVRMRZ
        mi_gold PFKFB3     89                    0                    0                    0                                                                                                                                                                                                                                                                                                            
        mi_gold  EPAS1     89                    0                    1                    1                                                                                                                                                                                                                                                                                              AHOUBRCZNHFOSL
        mi_gold   PDK4     89                    0                    0                    0                                                                                                                                                                                                                                                                                                            
