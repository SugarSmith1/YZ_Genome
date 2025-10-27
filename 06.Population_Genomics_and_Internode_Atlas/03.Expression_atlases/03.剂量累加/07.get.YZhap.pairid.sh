
# 唯一配对的基因对 9192
awk -F'\t' '($2 != "" && $3 != "" && $2 !~ /,|，/ && $3 !~ /,|，/)' cluster.YZhap_long.id> YZhap.pair.id