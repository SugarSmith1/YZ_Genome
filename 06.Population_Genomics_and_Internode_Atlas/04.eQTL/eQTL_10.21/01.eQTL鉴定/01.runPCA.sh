#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=hebhcnormal01

# Created on 2025-10-22
# @author: xiazhongqiang92@163.com

work_dir=/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/02.eQTL.10.22
cd ${work_dir}
mkdir 02.PCA
cd 02.PCA

~/miniconda3/bin/python3 ../pca_analysis_common.py \
    --input1 ${work_dir}/../01.quantitative/03.subgenome_long_gene_expre/YZhap.stage1.filter.tsv \
    --input2 ${work_dir}/../01.quantitative/03.subgenome_long_gene_expre/YZhap.stage2.filter.tsv \
    --input3 ${work_dir}/../01.quantitative/03.subgenome_long_gene_expre/YZhap.stage3.filter.tsv \
    --input4 ${work_dir}/../01.quantitative/03.subgenome_long_gene_expre/YZhap.stage4.filter.tsv \
    --output_dir ${work_dir}/02.PCA

for i in 1 2 3 4
do
    head -n 1 combined_pca_results.tsv > Period${i}_pca_results.tsv
    grep Period${i} combined_pca_results.tsv >> Period${i}_pca_results.tsv
done