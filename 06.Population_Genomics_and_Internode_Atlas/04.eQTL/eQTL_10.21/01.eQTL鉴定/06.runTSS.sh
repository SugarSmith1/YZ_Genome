#!/bin/bash
#SBATCH --job-name=submit_eqtl_jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --partition=hebhcnormal01

BASE_DIR="/public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/04.eQTL"
PEER_DIR="${BASE_DIR}/04.phe"
BED_DIR="${BASE_DIR}/04.bed"
PYTHON="/public/home/agis_xiazq/miniconda3/bin/python3"
GENE_BED="${BED_DIR}/filtered_YZhap_gene.bed"

for peer in 5 10 15 20 25 30 35 40; do
    for stage in 1 2 3 4; do
        JOB_NAME="eqtl_s${stage}_p${peer}"
        OUTPUT_BED="${PEER_DIR}/stage-${stage}_residuals-${peer}.bed"
        SCRIPT_PATH="${PEER_DIR}/run_${JOB_NAME}.sh"

        cat > "${SCRIPT_PATH}" <<EOF
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --partition=hebhcnormal01
#SBATCH --output=${PEER_DIR}/${JOB_NAME}.out
#SBATCH --error=${PEER_DIR}/${JOB_NAME}.err

cd ${PEER_DIR}

${PYTHON} ${BASE_DIR}/gene_TSS.py --peer_residuals ${PEER_DIR}/stage-${stage}_residuals-${peer}.tsv --gene_bed ${GENE_BED} --out_file ${OUTPUT_BED}

bgzip ${OUTPUT_BED}
tabix -p bed ${OUTPUT_BED}.gz
EOF

        sbatch "${SCRIPT_PATH}"
    done
done
