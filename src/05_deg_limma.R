# 05_deg_limma.R —— limma 差异分析（train 内主分析 + internal_val 复现检查）
# 输入：data/processed/GSE66360_gene.csv、data/processed/GSE66360_meta.csv
# 输出：results/tables/deg_train.csv、deg_train_relaxed.csv、deg_summary.csv、deg_internal.csv、deg_replication.csv
# 运行：C:/R/R-4.6.1/bin/Rscript.exe src/05_deg_limma.R

suppressPackageStartupMessages({
  library(limma)
})

# 项目根
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
root <- dirname(dirname(normalizePath(script_path, mustWork = FALSE)))
tab_dir <- file.path(root, "results", "tables")
log_dir <- file.path(root, "logs")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# 读取（check.names=FALSE 防止基因名被改写成 1-Mar/X1-Sep 之类）
expr <- read.csv(file.path(root, "data/processed/GSE66360_gene.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8", row.names = 1)
meta <- read.csv(file.path(root, "data/processed/GSE66360_meta.csv"),
                 check.names = FALSE, fileEncoding = "UTF-8")
colnames(meta) <- sub("^\ufeff", "", colnames(meta))
if (!"sample_id" %in% colnames(meta)) colnames(meta)[1] <- "sample_id"

expr <- expr[, meta$sample_id, drop = FALSE]

run_limma <- function(mat, grp) {
  g <- factor(grp, levels = c("control", "case"))
  design <- model.matrix(~0 + g)
  colnames(design) <- levels(g)
  fit <- lmFit(mat, design)
  cm <- makeContrasts(case - control, levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, number = Inf, adjust.method = "BH", sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  rownames(tt) <- NULL
  tt
}

# --- train 主分析
train_idx <- meta$role == "train"
expr_train <- expr[, train_idx, drop = FALSE]
grp_train <- meta$group[train_idx]
tt_train <- run_limma(expr_train, grp_train)
tt_train$sig <- (abs(tt_train$logFC) > 1) & (tt_train$adj.P.Val < 0.05)

tt_relaxed <- tt_train
tt_relaxed$sig <- (abs(tt_relaxed$logFC) > 0.5) & (tt_relaxed$adj.P.Val < 0.05)

n_sample_train <- sum(train_idx)
n_case <- sum(grp_train == "case")
n_control <- sum(grp_train == "control")
n_up <- sum(tt_train$sig & tt_train$logFC > 1)
n_down <- sum(tt_train$sig & tt_train$logFC < -1)
n_relaxed <- sum(tt_relaxed$sig)

write.csv(tt_train, file.path(tab_dir, "deg_train.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(tt_relaxed, file.path(tab_dir, "deg_train_relaxed.csv"), row.names = FALSE, fileEncoding = "UTF-8")

summary_df <- data.frame(
  dataset = "GSE66360", role = "train",
  n_sample = n_sample_train, n_case = n_case, n_control = n_control,
  threshold_main = "|logFC|>1 & adj.P.Val<0.05",
  n_up_main = n_up, n_down_main = n_down, n_sig_main = n_up + n_down,
  threshold_relaxed = "|logFC|>0.5 & adj.P.Val<0.05",
  n_sig_relaxed = n_relaxed
)
write.csv(summary_df, file.path(tab_dir, "deg_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("== train ==\n")
cat("样本", n_sample_train, "(case", n_case, "/ control", n_control, ")\n")
cat("主阈值 DEG：上调", n_up, "下调", n_down, "合计", n_up + n_down, "\n")
cat("放宽阈值 DEG：", n_relaxed, "\n")

if ((n_up + n_down) < 100) {
  cat("[STOP] train 主阈值 DEG =", n_up + n_down, "< 100，停止，等待人工决定\n")
  quit(status = 2)
}

# --- internal_val 复现检查（仅报告，不参与选特征）
internal_idx <- meta$role == "internal_val"
expr_internal <- expr[, internal_idx, drop = FALSE]
grp_internal <- meta$group[internal_idx]
tt_internal <- run_limma(expr_internal, grp_internal)
tt_internal$sig <- (abs(tt_internal$logFC) > 1) & (tt_internal$adj.P.Val < 0.05)
write.csv(tt_internal, file.path(tab_dir, "deg_internal.csv"), row.names = FALSE, fileEncoding = "UTF-8")

sig_train_genes <- tt_train$gene[tt_train$sig]
sig_internal_genes <- tt_internal$gene[tt_internal$sig]
overlap <- intersect(sig_train_genes, sig_internal_genes)
overlap_rate <- length(overlap) / mean(c(length(sig_train_genes), length(sig_internal_genes)))
dir_tr <- setNames(tt_train$logFC, tt_train$gene)
dir_in <- setNames(tt_internal$logFC, tt_internal$gene)
concordance <- if (length(overlap) > 0) mean(sign(dir_tr[overlap]) == sign(dir_in[overlap])) else NA

rep_df <- data.frame(
  metric = c("n_sig_train", "n_sig_internal", "n_overlap", "overlap_rate", "sign_concordance"),
  value = c(length(sig_train_genes), length(sig_internal_genes), length(overlap),
            round(overlap_rate, 4), round(concordance, 4))
)
write.csv(rep_df, file.path(tab_dir, "deg_replication.csv"), row.names = FALSE, fileEncoding = "UTF-8")

cat("== internal_val 复现 ==\n")
cat("internal 样本", sum(internal_idx), "\n")
cat("internal 显著", length(sig_internal_genes), "\n")
cat("重叠", length(overlap), "重叠率", round(overlap_rate, 4), "方向一致率", round(concordance, 4), "\n")

# 基因名抽查
all_genes <- rownames(expr)
probe_names <- c("MARCH1", "SEPT9", "MARCHF1", "SEPTIN9")
found <- probe_names[probe_names %in% all_genes]
bad <- grep("^(X[0-9]|[0-9]+-)", all_genes, value = TRUE)
cat("基因名抽查：", paste(found, collapse = ", "), "\n")
cat("疑似被改写基因名数：", length(bad), "\n")

# top10 上调/下调（写入日志）
top_up_tab <- tt_train[tt_train$sig & tt_train$logFC > 1, ]
top_up <- top_up_tab[order(-top_up_tab$logFC), , drop = FALSE]
top_up <- top_up$gene[1:min(10, nrow(top_up))]
top_down_tab <- tt_train[tt_train$sig & tt_train$logFC < -1, ]
top_down <- top_down_tab[order(top_down_tab$logFC), , drop = FALSE]
top_down <- top_down$gene[1:min(10, nrow(top_down))]

log_lines <- c(
  "# P3 limma 差异分析日志", "",
  paste0("- limma 版本：", as.character(packageVersion("limma"))),
  paste0("- train 样本数：", n_sample_train, "（case ", n_case, " / control ", n_control, "）"),
  paste0("- 主阈值：|logFC|>1 & adj.P.Val<0.05 -> 上调 ", n_up, " / 下调 ", n_down, " / 合计 ", n_up + n_down),
  paste0("- 放宽阈值：|logFC|>0.5 & adj.P.Val<0.05 -> ", n_relaxed),
  paste0("- internal_val 样本数：", sum(internal_idx)),
  paste0("- internal 显著：", length(sig_internal_genes)),
  paste0("- train/internal 重叠：", length(overlap), "，重叠率 ", round(overlap_rate, 4), "，方向一致率 ", round(concordance, 4)),
  paste0("- 基因名抽查：", paste(found, collapse = ", "), "；疑似被改写数：", length(bad)),
  paste0("- top10 上调：", paste(top_up, collapse = ", ")),
  paste0("- top10 下调：", paste(top_down, collapse = ", ")),
  ""
)
writeLines(log_lines, file.path(log_dir, "deg_log.md"), useBytes = TRUE)

cat("[05_deg_limma] 完成\n")
