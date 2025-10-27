# -*- coding: utf-8 -*-

import pandas as pd
import click

@click.command()
@click.option("--peer_residuals", type=click.Path(exists=True), required=True, help="PEER矫正后的残差数据文件路径（CSV格式）。")
@click.option("--gene_bed", type=click.Path(exists=True), required=True, help="基因的BED文件路径（包含基因坐标信息）。")
@click.option("--out_file", type=click.Path(), required=True, help="输出文件路径（CSV格式）。")
def main(peer_residuals, gene_bed, out_file):
    """
    将PEER矫正后的数据处理为tensorQTL输入格式。
    """
    # 读取PEER矫正后的残差数据
    peer_residuals_df = pd.read_csv(peer_residuals, header=0, index_col=0, sep="\t")

    # 读取基因BED文件
    gene_bed_df = pd.read_csv(gene_bed, index_col=4, sep="\t", header=None)
    gene_bed_df[0] = gene_bed_df[0].apply(lambda x: int(x))  
    gene_bed_df['TSS'] = gene_bed_df.apply(lambda x: x[1] if x[3] == "+" else x[2], axis=1)  # 计算转录起始位点（TSS）

    # 处理数据并生成tensorQTL输入格式
    gene_info = []
    expression_info = []
    for gene_id, express_data in peer_residuals_df.T.iterrows():
        if gene_id not in gene_bed_df.index:
            # 跳过scaffold上的基因
            continue
        else:
            chrom, start, end, strand, tss = gene_bed_df.loc[gene_id]
            gene_info.append([int(chrom), tss, tss + 1, gene_id])  # 添加基因坐标信息
            expression_info.append(express_data.values)  # 添加表达数据

    # 将基因信息和表达数据合并为DataFrame
    gene_info_df = pd.DataFrame(gene_info)
    expression_info_df = pd.DataFrame(expression_info)
    out_df = pd.concat([gene_info_df, expression_info_df], axis=1)
    out_df.columns = ['#chr', 'start', 'end', 'phenotype'] + list(peer_residuals_df.T.columns)

    # 按照染色体和起始位置排序
    out_df = out_df.sort_values(by=['#chr', 'start'])

    # 保存结果
    out_df.to_csv(out_file, header=True, index=False, sep="\t")
    print(f"tensorQTL输入格式的数据已保存至: {out_file}")

if __name__ == "__main__":
    main()