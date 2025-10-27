#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=hebhcnormal01


echo started at `date` @`hostname`
module load apps/samtools/1.9/gcc-7.3.1


export LD_LIBRARY_PATH=/public/software/apps/openssl/1.1.1l/lib/:$LD_LIBRARY_PATH

~/software/samtools-1.20/bin/samtools merge -o map_hifi.bam map_hifi_1.bam map_hifi_2.bam map_hifi_3.bam -@ 64
~/software/samtools-1.20/bin/samtools sort map_hifi.bam -o map_hifi.sorted.bam -@ 64


echo ended at `date`