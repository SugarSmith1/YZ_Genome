#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript script.R <input_expr_file> <min_samples_expr> <output_file>\n")
}

input_file <- args[1]
min_samples_expr <- as.integer(args[2])
output_file <- args[3]

if (!file.exists(input_file)) {
  stop(paste("Input expression file does not exist:", input_file))
}

# 读取表达矩阵，行为基因，列为样本
expr <- read.delim(input_file, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cat("Read expression matrix:", dim(expr), "\n")

# 定义表达阈值 TPM > 0.01 为表达
expr_binary <- expr > 0.01

# 统计每个基因在多少样本中表达
expressed_counts <- rowSums(expr_binary)

# 过滤基因：保留在至少 min_samples_expr 样本中表达的基因
expr_filtered <- expr[expressed_counts >= min_samples_expr, ]
cat("Genes retained after filtering:", nrow(expr_filtered), "\n")

if (nrow(expr_filtered) == 0) {
  stop("No genes pass the expression filter. Please lower the threshold or check data.")
}

# 保存过滤后的原始表达矩阵
write.table(expr_filtered, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("Filtered expression matrix saved to:", output_file, "\n")
