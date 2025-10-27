#!/usr/bin/env Rscript

library(ggplot2)

md <- read.table("hifi_depth.chr.txt", header=F, sep="\t")

names(md) <- c("chrom", "chromStart", "chromEnd", "readCount", "meanCoverage", "sampleName")

md$chrom <- factor(md$chrom, levels=c(
    "chr1_1", "chr1_2", "chr1_3", "chr1_4", "chr1_5", "chr1_6", "chr1_7", "chr1_8", "chr1_9", "chr1_10", "chr1_11", "chr1_12",
    "chr2_1", "chr2_2", "chr2_3", "chr2_4", "chr2_5", "chr2_6", "chr2_7", "chr2_8", "chr2_9", "chr2_10", "chr2_11", "chr2_12",
    "chr3_1", "chr3_2", "chr3_3", "chr3_4", "chr3_5", "chr3_6", "chr3_7", "chr3_8", "chr3_9", "chr3_10", "chr3_11", "chr3_12",
    "chr4_1", "chr4_2", "chr4_3", "chr4_4", "chr4_5", "chr4_6", "chr4_7", "chr4_8", "chr4_9", "chr4_10", "chr4_11", "chr4_12",
    "chr5_1", "chr5_2", "chr5_3", "chr5_4", "chr5_5", "chr5_6", "chr5_7", "chr5_8", "chr5_9",
    "chr6_1", "chr6_2", "chr6_3", "chr6_4", "chr6_5", "chr6_6", "chr6_7", "chr6_8", "chr6_9", "chr6_10", "chr6_11", "chr6_12",
    "chr7_1", "chr7_2", "chr7_3", "chr7_4", "chr7_5", "chr7_6", "chr7_7", "chr7_8", "chr7_9", "chr7_10", "chr7_11", "chr7_12",
    "chr8_1", "chr8_2", "chr8_3", "chr8_4", "chr8_5", "chr8_6", "chr8_7", "chr8_8", "chr8_9",
    "chr9_1", "chr9_2", "chr9_3", "chr9_4", "chr9_5", "chr9_6", "chr9_7", "chr9_8", "chr9_9", "chr9_10", "chr9_11", "chr9_12",
    "chr10_1", "chr10_2", "chr10_3", "chr10_4", "chr10_5", "chr10_6", "chr10_7", "chr10_8", "chr10_9", "chr10_10", "chr10_11"
    ))

pdf("depth.pdf", width=10, height=18)

ggplot(md, aes(x=(chromEnd+chromStart)/2000000, y=meanCoverage)) + geom_line() + theme_test() + ylim(0, 70) + 
    facet_wrap(~chrom, scale="free_x", ncol=6) + xlab("Position (Mb)") + ylab("Sequencing depth")

dev.off()