work_dir=/public/home/agis_xiazq/project/01.YZ081609/05.Population_Profiling/08.bias-eQTL
cd ${work_dir}
mkdir 02.PCA
cd 02.PCA

python3 ../pca_analysis_common.py \
    --input1 ${work_dir}/01.exp/updated_so_Hap_genes_stage1.tsv \
    --input2 ${work_dir}/01.exp/updated_so_Hap_genes_stage2.tsv \
    --input3 ${work_dir}/01.exp/updated_so_Hap_genes_stage3.tsv \
    --input4 ${work_dir}/01.exp/updated_so_Hap_genes_stage4.tsv \
    --output_dir ${work_dir}/02.PCA

for i in 1 2 3 4
do
    head -n 1 combined_pca_results.tsv > Period${i}_pca_results.tsv
    grep Period${i} combined_pca_results.tsv >> Period${i}_pca_results.tsv
done