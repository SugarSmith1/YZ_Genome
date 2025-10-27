#!/bin/bash

exp_data=/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/01.quantitative/03.subgenome_long_gene_expre
work_dir=/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/02.eQTL.10.22
cd ${work_dir}
mkdir 03.peer_interface
cd 03.peer_interface
# 为每个stage创建一个独立的sbatch任务
for j in 5 10 15 20 25 30 35 40
do
for i in 1 2 3 4
do
    job_script="peer_stage_${i}_${j}.sh"

    # 生成sbatch脚本文件
    cat > $job_script <<EOF
#!/bin/bash
#SBATCH --job-name=peer_stage_${i}_${j}
#SBATCH --output=logs/peer_stage_${i}_${j}.out
#SBATCH --error=logs/peer_stage_${i}_${j}.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --partition=hebhcnormal01

source ~/miniconda3/bin/activate r_peer
cd ${work_dir}/03.peer_interface
/public/home/agis_xiazq/miniconda3/envs/r_peer/bin/Rscript /public/home/agis_xiazq/Script/eQTL/peer.r ${exp_data}/YZhap.stage${i}.filter.tsv ${j} results/stage_${i}
EOF

    # 创建日志目录（如果不存在）
    mkdir -p logs

    # 提交生成的sbatch脚本
    sbatch $job_script
    sleep 5
    echo "Submitted stage ${i}_${j} to sbatch."
done
done