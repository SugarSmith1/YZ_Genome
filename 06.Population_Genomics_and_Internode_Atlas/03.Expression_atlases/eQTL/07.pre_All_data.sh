

# 准备所有eQTL数据，并迁移至小服务器

GWAS_SOURCE=~/project/000.Nature/13.GWAS/01.gwas_vcf-2/final-gwas
EXPR_SOURCE=/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/04.eQTL/04.phe
PCA_SOURCE=/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/02.pca

mkdir 07.pre.all.data && cd 07.pre.all.data
cp ${GWAS_SOURCE}/GWAS* .
cp ${EXPR_SOURCE}/*bed.gz* .

for i in 1 2 3 4
do
    grep "Cul" ${PCA_SOURCE}/Period${i}_pca_results.tsv | \
    awk '{print $1"\t"$1"\t"$2"\t"$3"\t"$4}' > PCA_qcovar.Stage${i}.txt
done
