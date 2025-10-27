#!/bin/bash
# QTL分析调用脚本

source ~/miniconda3/bin/activate tensorqtl_env

# 定义参数
stages=(1 2 3 4)
modes=("p" "t")  # p=cis-eQTL, t=trans-eQTL
factors=(5 10 15 20 25 30 35 40)

# 创建输出目录（如果不存在）
output_base="/data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/05.pair_eqtl"
mkdir -p $output_base

# 遍历所有组合
for stage in "${stages[@]}"; do
  for mode in "${modes[@]}"; do
    for factor in "${factors[@]}"; do
      
      # 构建输出文件前缀
      output_prefix="${output_base}/${stage}_${mode}_${factor}"
      
      echo "=========================================="
      echo "Processing: Stage=${stage}, Mode=${mode}, Factor=${factor}"
      echo "Output: ${output_prefix}_*"
      echo "=========================================="
      
      # 运行Python脚本
      python3 qtl_analysis.py \
        --expression_bed /data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data/07.pre.all.data/stage-${stage}_residuals-${factor}.bed.gz \
        --covariates_file /data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data/07.pre.all.data/PCA_qcovar.Stage${stage}.txt \
        --outfile $output_prefix \
        --mode $mode
      
      # 检查输出文件
      if [ "$mode" == "p" ] || [ "$mode" == "both" ]; then
        if [ -f "${output_prefix}_cis_lead_snps.txt" ]; then
          line_count=$(wc -l < "${output_prefix}_cis_lead_snps.txt")
          echo "Cis-eQTL output: ${output_prefix}_cis_lead_snps.txt ($((line_count-1)) significant genes)"
        fi
      fi
      
      if [ "$mode" == "t" ] || [ "$mode" == "both" ]; then
        if [ -f "${output_prefix}_trans_all_significant.txt" ]; then
          line_count=$(wc -l < "${output_prefix}_trans_all_significant.txt")
          echo "Trans-eQTL output: ${output_prefix}_trans_all_significant.txt ($((line_count-1)) significant pairs)"
        fi
      fi
      
      echo ""
      
    done
  done
done

echo "所有QTL分析完成！"