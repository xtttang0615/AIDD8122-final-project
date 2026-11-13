# 输入文件：无
# 输出文件：logs/env_report.md（追加）
# 运行命令：Rscript src/00_check_env.R
# 说明：用 requireNamespace 逐个检查 Bioconductor 包，输出 TRUE/FALSE 并追加到报告。

# 从 --file= 参数推导脚本路径，进而得到项目根目录
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")
} else {
  script_path <- normalizePath(".", winslash = "/")
}
root <- dirname(dirname(script_path))

pkg_names <- c("limma", "GEOquery", "AnnotationDbi", "hgu133plus2.db", "Biobase")

res <- data.frame(package = pkg_names, available = NA, version = NA, stringsAsFactors = FALSE)
for (i in seq_along(pkg_names)) {
  p <- pkg_names[i]
  ok <- requireNamespace(p, quietly = TRUE)
  res$available[i] <- ok
  res$version[i] <- if (ok) as.character(packageVersion(p)) else "MISSING"
}

report_path <- file.path(root, "logs", "env_report.md")

out <- c(
  "",
  "## 4. R packages (Bioconductor)",
  "",
  "| package | available | version |",
  "|---|---|---|"
)
for (i in seq_len(nrow(res))) {
  out <- c(out, sprintf("| %s | %s | %s |",
                        res$package[i],
                        toupper(as.character(res$available[i])),
                        res$version[i]))
}

# 追加写入报告
con <- file(report_path, open = "at", encoding = "UTF-8")
writeLines(out, con = con)
close(con)

# 打印到控制台
cat(paste(out, collapse = "\n"), "\n")

# 缺包记录到 setup/missing_packages.txt
setup_dir <- file.path(root, "setup")
if (!dir.exists(setup_dir)) dir.create(setup_dir, recursive = TRUE)
missing_txt <- file.path(setup_dir, "missing_packages.txt")
r_missing <- res$package[!res$available]
if (length(r_missing) > 0) {
  line <- paste("R packages missing:", paste(r_missing, collapse = ", "))
} else {
  line <- "R packages missing: none"
}
cat(line, "\n", file = missing_txt, append = TRUE)