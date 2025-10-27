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

# 定义自定义分位数归一化函数
quantile_normalisation <- function(df){
  df_rank <- apply(df, 2, rank, ties.method="min")
  df_sorted <- apply(df, 2, sort)
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean) my_mean[my_index]
  df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}

# 定义表达阈值TPM > 0.01为表达
expr_binary <- expr > 0.01

# 统计每个基因在多少样本中表达
expressed_counts <- rowSums(expr_binary)

# 过滤基因：保留在至少 min_samples_expr 样本中表达的基因
expr_filtered <- expr[expressed_counts >= min_samples_expr, ]
cat("Genes retained after filtering:", nrow(expr_filtered), "\n")

if (nrow(expr_filtered) == 0) {
  stop("No genes pass the expression filter. Please lower the threshold or check data.")
}

# log2(x + 1)转换
expr_log <- log2(expr_filtered + 1)

# 使用自定义函数进行分位数归一化
expr_norm <- quantile_normalisation(expr_log)

cat("Normalized expression matrix dimension:", dim(expr_norm), "\n")

# 保存结果，带行名和列名
# write.table(expr_norm, file = output_file, sep = "\t", quote = FALSE, col.names = NA)
write.table(expr_norm, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

cat("Normalized data saved to:", output_file, "\n")