# 加载 PEER 库
library(peer)

# 从命令行参数获取文件路径和因子数量
args <- commandArgs(trailingOnly = TRUE)

# 参数解析
if (length(args) < 3) {
  stop("请提供以下参数：\n1. 表达数据文件路径\n2. 因子数量\n3. 输出目录路径")
}

expr_file <- args[1]         # 表达数据文件路径
num_factors <- as.integer(args[2])  # 因子数量
output_dir <- args[3]        # 输出目录路径

# 读取表达数据
if (!file.exists(expr_file)) {
  stop(paste("表达数据文件不存在：", expr_file))
}
expr <- read.delim(expr_file, header = TRUE, row.names = 1, sep = "\t")
cat("输入数据的维度：", dim(expr), "\n")

# 创建 PEER 模型
model <- PEER()

# 设置表达数据
PEER_setPhenoMean(model, as.matrix(t(expr)))

# 设置因子数量
PEER_setNk(model, num_factors)
cat("设置的因子数量：", PEER_getNk(model), "\n")

# 更新模型
PEER_update(model)

# 获取因子和残差
factors <- as.data.frame(t(PEER_getX(model)))
residuals <- as.data.frame(t(PEER_getResiduals(model)))

# 设置列名
colnames(factors) <- colnames(expr)
colnames(residuals) <- colnames(expr)

# 创建输出目录
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 保存因子和残差
factors_file <- file.path(output_dir, paste0('factors_', num_factors, '.txt'))
residuals_file <- file.path(output_dir, paste0('residuals_', num_factors, '.txt'))

write.table(factors, file = factors_file, sep = '\t', quote = FALSE, col.names = FALSE)
write.table(residuals, file = residuals_file, sep = '\t', quote = FALSE, col.names = FALSE)

cat("PEER 模型分析完成，结果已保存至：\n",
    "- 因子文件：", factors_file, "\n",
    "- 残差文件：", residuals_file, "\n")