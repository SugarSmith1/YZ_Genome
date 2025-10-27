#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=hebhcnormal01

work_dir=/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/02.eQTL.10.22
cd ${work_dir}
mkdir 04.phe
cd 04.phe

for peer in 5
do
        for stage in 1 2 3 4
        do
        ~/miniconda3/bin/python3 ../peer_RINT.py --peer_file /public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/02.eQTL.10.22/03.peer_interface/results/stage_${stage}/residuals_${peer}.txt --expre_file /public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/01.quantitative/03.subgenome_long_gene_expre/YZhap.stage${stage}.filter.tsv --output ${work_dir}/04.phe/stage-${stage}_residuals-${peer}.tsv
        done
done