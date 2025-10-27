
for i in 1 2 3 4
do
python3 extract_hap_gene_expression.py --cluster-file ~/project/02.YZ/R2.SubGenome/04.allele_defined_subgenome/03.get-allele-table/long_id/cluster.YZhap_long.id --so-expr-file /public/home/agis_xiazq/project/02.YZ/R4.PopVar/quantitative/stage${i}_sum_expression_data_so_sums.tsv --ss-expr-file /public/home/agis_xiazq/project/02.YZ/R4.PopVar/quantitative/stage${i}_sum_expression_data_ss_sums.tsv --output-file YZhap.stage${i}.tsv
done