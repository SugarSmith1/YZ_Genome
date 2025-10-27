#!/usr/bin/env bash
# count_by_stage.sh  — four stages × (scaf/SO/SS)



lists=(scaf SO SS)      # 需已有：YZ.scaf.gene.list  YZ.SO.gene.list  YZ.SS.gene.list
echo -e "Stage\tType\tTotal\tExpressed\tFiltered" > expression_summary.tsv

for s in 1 2 3 4; do
  orig="stage${s}_ori_expression_data.tsv"
  filt="stage${s}.filter.tsv"

  cut -f1 "$orig" | tail -n +2 > all.tmp
  cut -f1 "$filt" | tail -n +2 > filt.tmp
  grep -Fvxf filt.tmp all.tmp > expr.tmp

  for t in "${lists[@]}"; do
    list="YZ.${t}.gene.list"
    total=$(wc -l < "$list")
    expr=$(grep -Fxf "$list" expr.tmp | wc -l)
    filn=$(grep -Fxf "$list" filt.tmp | wc -l)
    echo -e "stage${s}\t${t}\t${total}\t${expr}\t${filn}" >> expression_summary.tsv
  done

  rm -f all.tmp filt.tmp expr.tmp
done

cat expression_summary.tsv
