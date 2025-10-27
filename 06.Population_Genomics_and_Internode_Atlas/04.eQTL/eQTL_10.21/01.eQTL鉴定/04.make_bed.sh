
## 获取总的YZhap，gene-bed
awk '{
  if ($3 == "gene") {
    split($9, a, ";");
    split(a[1], b, "=");
    print $1"\t"$4"\t"$5"\t"$7"\t"b[2]
  }
}' /public/home/agis_xiazq/project/02.YZ/R2.SubGenome/02.genome_annotation/YZhap.Chrom.gff3 > YZhap_gene.bed

