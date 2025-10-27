
## 获取总的YZhap，gene-bed
hap_dir=/public/home/agis_xiazq/project/000.Nature/01.genome_assembly/02.chrom/01.YZ081609/YZhap
awk '{
  if ($3 == "gene") {
    split($9, a, ";");
    split(a[1], b, "=");
    print $1"\t"$4"\t"$5"\t"$7"\t"b[2]
  }
}' /public/home/agis_xiazq/project/000.Nature/01.genome_assembly/02.chrom/01.YZ081609/YZhap/YZhap.Chrom.gff3 > YZhap_gene.bed

## 提取pair_gene

awk 'NR==FNR{genes[$1]; genes[$2]; next} $5 in genes' <(awk '{print $2,$3}' /public/home/agis_xiazq/project/02.YZ/R4.PopVar/pop_expression/01.quantitative/02.subgenome_gene_expre/YZhap.pair.id | tr ' ' '\n') YZhap_gene.bed > filtered_YZhap_gene.bed

