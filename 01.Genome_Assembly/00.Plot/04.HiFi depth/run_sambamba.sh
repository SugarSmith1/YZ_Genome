#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --partition=hebhcnormal01


echo started at `date` @`hostname`

export LD_LIBRARY_PATH=/public/software/apps/openssl/1.1.1l/lib/:$LD_LIBRARY_PATH

# sambamba index map_hifi.sorted.bam -t 32
sambamba depth window map_hifi.sorted.bam -w 1000000 -F 'mapping_quality >= 0 and not duplicate and not failed_quality_control' -t 32 -o hifi_depth.txt


echo ended at `date`