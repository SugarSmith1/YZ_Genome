## 提取同源基因对的表达
for i in 1 2 3 4
do
python3 extract_hap_gene_expression.py --cluster-file YZhap.pair.id --so-expr-file /public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/01.quantitative/02.subgenome_gene_expre/stage${i}_sum_expression_data_so_sums.tsv --ss-expr-file /public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/01.quantitative/02.subgenome_gene_expre/stage${i}_sum_expression_data_ss_sums.tsv --output-file YZhap_pair.stage${i}.tsv
done