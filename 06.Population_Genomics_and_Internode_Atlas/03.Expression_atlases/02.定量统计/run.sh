
# 对YZ所有基因进行过滤
bash run.filter.sh

# 准备亚基因组拆分的基因集
gene_list=~/project/000.Nature/03.lineage_decomposition/01.YZ/04.lineage_result/YZ081609-noStep/YZ-SO-SS-gene-debug3.txt
grep scaffold ${gene_list} |cut -f1 > YZ.scaf.gene.list
grep SO ${gene_list} |cut -f1 > YZ.SO.gene.list
grep SS ${gene_list} |cut -f1 > YZ.SS.gene.list

# 统计过滤和非过滤的非类基因占比


