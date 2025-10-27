#!/bin/bash
# 小服务器

source ~/miniconda3/bin/activate tensorqtl_env
# 定义阶段、模式和因子列表
stages=(1 2 3 4)
#modes=("p" "t" "n")
#modes=("p" "t" "n")
modes=("p" "t")
factors=(5 10 15 20 25 30 35 40)
# factors=(5 10 40)
# 遍历每个阶段、模式和因子数量的组合
for stage in "${stages[@]}"; do
  for mode in "${modes[@]}"; do
    for factor in "${factors[@]}"; do
      python3 QTL_mapping.py \
        --expression_bed /data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data/07.pre.all.data/stage-${stage}_residuals-${factor}.bed.gz \
        --covariates_file /data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data/07.pre.all.data/PCA_qcovar.Stage${stage}.txt \
        --outfile /data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/05.pair_eqtl/${stage}_${mode}_${factor} \
        --mode $mode
      #python3 mv_parquet_txt.py ./
      
      # 打印运行信息
      echo "Processed: Stage=${stage}, Mode=${mode}, Factor=${factor}"
    done
  done
done

echo "所有任务已完成！"