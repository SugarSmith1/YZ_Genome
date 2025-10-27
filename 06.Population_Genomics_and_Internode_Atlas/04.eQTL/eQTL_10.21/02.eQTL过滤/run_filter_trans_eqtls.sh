python3 filter_trans_eqtls.py \
  --input stage1_t_5_trans_all_significant.txt \
  --output stage1_t_5_trans_strict.txt \
  --qval-threshold 0.01 \
  --pval-threshold 1e-10 \
  --effect-size-threshold 0.2 \
  --min-maf 0.1 \
  --snps-per-gene 3 \
  --filter-hotspots \
  --filter-cis-acting

# >5mb 或者不同染色体
python filter_true_trans_eqtls.py \
  --input stage1_t_5_trans_all_significant.txt \
  --output stage1_t_5_true_trans.txt \
  --gene-bed /path/to/gene_positions.bed \
  --distance-threshold 5000000 \
  --qval-threshold 0.05 \
  --pval-threshold 1e-8

