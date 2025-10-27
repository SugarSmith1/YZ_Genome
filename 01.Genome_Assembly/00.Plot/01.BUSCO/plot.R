#!/usr/bin/env Rscript

library(ggplot2)

md <- read.table("busco_count.txt", header=F, sep=" ")
names(md) <- c("Count", "Ploidy", "Version")
md$Version <- factor(md$Version, levels=c("final", "Nature", "NC", "NG"))


pdf("busco.pdf", width=6, height=4)
ggplot(md, aes(x=Ploidy, y=Count)) + geom_line(aes(color=Version)) + geom_point(aes(color=Version)) + 
    theme_test() + scale_color_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5388FF"))
dev.off()