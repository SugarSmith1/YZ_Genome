#!/usr/bin/env Rscript

library(ggplot2)
library(patchwork)
library(ggpubr)

md <- read.table("/Users/xiazhongqiang/博士后/02.项目/01.甘蔗高蔗糖的多基因组学与遗传学/02.云蔗手稿/01.YZ081609_Code/01.文章绘图/绘制组装优势性/results.txt", header=F, sep="\t")

names(md) <- c("Assembly", "Query", "Query_len", "Reference", "Reference_len", "Orientation", "Chained_aln_len", "Chained_aln_num", "Percent_identity", "Gap_compressed_percent_identity")

md$Assembly <- factor(md$Assembly, levels=c("YZ", "Nature", "NC", "NG"))

md$Reference <- factor(md$Reference, levels=c(
        "GWHBRAN00000001", "GWHBRAN00000002", "GWHBRAN00000003", "GWHBRAN00000004", "GWHBRAN00000005",
        "GWHBRAN00000006", "GWHBRAN00000007", "GWHBRAN00000008", "GWHBRAN00000009", "GWHBRAN00000010"))

comparisons <- list(c("YZ", "Nature"), c("YZ", "NC"), c("YZ", "NG"))

# 修改1：调整 p1 的 y 轴标签，仅左侧子图显示
p1 <- ggplot(md, aes(x=Assembly, y=Chained_aln_len/1000000)) + 
    geom_boxplot(aes(fill=Assembly)) + 
    facet_wrap(~Reference, ncol=5, strip.position = "top") +  # 调整分面位置
    stat_compare_means(comparisons=comparisons, method="wilcox.test", paired=FALSE, label="p.format") + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    ylab("Total chained alignment length (Mb)") + 
    ggtitle("Total chained alignment length") +
    theme(
        axis.title.y = element_text(),  # 保留 y 轴标题
        axis.text.y = element_text(),   # 保留 y 轴刻度
        strip.background = element_blank(),  # 可选：移除分面标签背景
        panel.spacing = unit(0.5, "lines")   # 调整子图间距
    )

# 修改2：调整 p2 的 y 轴标签，仅左侧子图显示
p2 <- ggplot(md, aes(x=Assembly, y=Gap_compressed_percent_identity)) + 
    geom_boxplot(aes(fill=Assembly)) + 
    facet_wrap(~Reference, ncol=5, strip.position = "top") +  # 调整分面位置
    stat_compare_means(comparisons=comparisons, method="wilcox.test", paired=FALSE, label="p.format") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    ylab("Gap compressed percent identity") + 
    ggtitle("Gap compressed percent identity") +
    theme(
        axis.title.y = element_text(),  # 保留 y 轴标题
        axis.text.y = element_text(),   # 保留 y 轴刻度
        strip.background = element_blank(),  # 可选：移除分面标签背景
        panel.spacing = unit(0.5, "lines")   # 调整子图间距
    )

# 使用 patchwork 组合图形时，进一步调整 y 轴标签显示
pdf("results.pdf", width=10, height=10)

patchwork <- p1 + p2 + plot_layout(guides='collect', ncol=1) + plot_annotation(tag_levels="a")

# 移除右侧子图的 y 轴标签
final_plot <- patchwork & 
    theme_test() & 
    theme(
        legend.position = "bottom",
        legend.margin = margin(c(0, 0, 0, -2), unit = "cm"),
        axis.title.y.right = element_blank(),  # 移除右侧 y 轴标题
        axis.text.y.right = element_blank(),   # 移除右侧 y 轴刻度
        axis.ticks.y.right = element_blank()   # 移除右侧 y 轴刻度线
    ) &
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5388FF"))

print(final_plot)

dev.off()
