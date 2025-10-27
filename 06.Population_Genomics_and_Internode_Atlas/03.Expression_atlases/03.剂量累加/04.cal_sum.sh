

# 根据等位基因表，计算每个簇的总表达情况并区分亚基因组

for i in 1 2 3 4
do
~/miniconda3/bin/python3 cal_sum_uniq.py --cluster-file /public/home/agis_xiazq/project/02.YZ/R2.SubGenome/03.allele_defined_subgenome/03.get-allele-table/genes_output.tsv --expression-file /public/home/agis_xiazq/project/000.Nature/12-1.transcriptome/09.YZ_genome_kallisto_pop/02.gene.expre/stage${i}_ori_expression_data.tsv --output-prefix stage${i}_sum_expression_data
done
