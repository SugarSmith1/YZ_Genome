
## 过滤基因
# TPM > 0.01 定义为表达，保留至少50个样本表达的基因；
for i in 1 2 3 4
do
Rscript prepare_gene_expression_onlyfilter.R stage${i}_ori_expression_data.tsv 50 stage${i}.filter.tsv
done
