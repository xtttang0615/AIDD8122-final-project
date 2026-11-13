# 05b_limma_fold.R —— 供 07_signature.py 调用：在单个外层折的训练部分内跑 limma
# 用法：Rscript 05b_limma_fold.R <expr_csv> <labels_csv> <fold_no> <out_csv>
#   expr_csv   ：训练集表达矩阵（行=基因，列=GSM，UTF-8）
#   labels_csv ：列含 sample_id, group, fold（fold=该样本所属留出折号 1..5）
#   fold_no    ：当前外层折号，训练部分 = fold != fold_no
#   out_csv    ：输出 limma 全表（gene/logFC/AveExpr/t/P.Value/adj.P.Val/B）

suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
expr_csv <- args[1]
label_csv <- args[2]
fold_no <- as.integer(args[3])
out_csv <- args[4]

expr <- read.csv(expr_csv, check.names = FALSE, fileEncoding = "UTF-8", row.names = 1)
lab <- read.csv(label_csv, check.names = FALSE, fileEncoding = "UTF-8")

sub <- lab[lab$fold != fold_no, , drop = FALSE]
expr <- expr[, sub$sample_id, drop = FALSE]

g <- factor(sub$group, levels = c("control", "case"))
design <- model.matrix(~0 + g)
colnames(design) <- levels(g)

fit <- lmFit(expr, design)
cm <- makeContrasts(case - control, levels = design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)
tt <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")
tt$gene <- rownames(tt)
tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]

write.csv(tt, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("fold", fold_no, ": n_train =", ncol(expr), ", n_genes =", nrow(expr), "\n")
