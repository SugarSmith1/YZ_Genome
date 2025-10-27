## 过滤基因
# TPM > 0.01 定义为表达，保留至少30个样本表达的基因；对过滤后的数据进行 log2(x+1) 转换；进行分位数归一化
for i in 1 2 3 4
do
Rscript prepare_gene_expression.R YZhap.stage${i}.tsv 50 YZhap.stage${i}.filter.tsv
done

## 过滤基因
# TPM > 0.01 定义为表达，保留至少30个样本表达的基因；对过滤后的数据进行 log2(x+1) 转换；进行分位数归一化
for i in 1 2 3 4
do
Rscript prepare_gene_expression.R YZhap_pair.stage${i}.tsv 50 YZhap_pair.stage${i}.filter.tsv
done