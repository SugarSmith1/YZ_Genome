# -*- coding: utf-8 -*-
'''
Created on 2024-12-14 16:06:20

@author: xiazhongqiang92@163.com
'''

import os
import sys
import torch
import click
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
from statsmodels.stats.multitest import multipletests

# 设置 CUDA 设备
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Using device: {device}')

# 硬编码 PLINK 文件路径
PLINK_PREFIX_PATH = "/data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data/07.pre.all.data/GWAS"

def load_expression_data(expression_bed):
    """加载表达量数据."""
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    print(f"Loaded expression data: {phenotype_df.shape}")
    return phenotype_df, phenotype_pos_df

def load_covariates(covariates_file, samples):
    # 加载协变量文件
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0, header=None)
    # 删除不需要的列（例如第 1 列）
    covariates_df = covariates_df.drop([1], axis=1)
    # 设置索引名称和列名
    covariates_df.index.name = 'id'
    covariates_df.columns = ['PC1', 'PC2', 'PC3']
    # 筛选匹配样本
    covariates_df = covariates_df.loc[samples]
    return covariates_df


def perform_cis_analysis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, outfile, nperm, maf_threshold, window):
    cis_df = cis.map_cis(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
        nperm=nperm, maf_threshold=maf_threshold, window=window, seed=2022
    )
    # 替代calculate_qvalues的FDR校正
    pvals = cis_df['pval_nominal'].values
    _, qvals, _, _ = multipletests(pvals, method='fdr_bh')
    cis_df['qval'] = qvals
    
    cis_df.to_csv(outfile, header=True, index=True, sep="\t")
    print(f"Cis-eQTL analysis results saved to {outfile}")

def perform_nominal_mapping(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, outfile, maf_threshold, window):
    """执行nominal映射分析."""
    cis.map_nominal(
        genotype_df, variant_df, phenotype_df, phenotype_pos_df,
        prefix=outfile, covariates_df=covariates_df,
        maf_threshold=maf_threshold, window=window, output_dir='.',
        write_top=True, write_stats=True
    )
    print(f"Nominal mapping results saved with prefix {outfile}")

def perform_trans_analysis(genotype_df, phenotype_df, covariates_df, outfile, pval_threshold, maf_threshold):
    """执行trans-QTL分析."""
    trans_df = trans.map_trans(
        genotype_df, phenotype_df, covariates_df,
        return_sparse=True, pval_threshold=pval_threshold, maf_threshold=maf_threshold, batch_size=20000
    )
    trans_df.to_csv(outfile, header=True, index=False)
    print(f"Trans-QTL analysis results saved to {outfile}")

@click.command()
@click.option('--expression_bed', type=click.Path(exists=True), required=True, help="Path to expression BED file.")
@click.option('--covariates_file', type=click.Path(exists=True), required=True, help="Path to covariates file.")
@click.option('--outfile', type=click.Path(), required=True, help="Path to save the output results.")
@click.option('--mode', type=click.Choice(['p', 'n', 't']), required=True, help="Mode of operation: 'p' for cis-eQTL, 'n' for nominal mapping, 't' for trans-QTL mapping.")
@click.option('--nperm', type=int, default=1000, show_default=True, help="Number of permutations for cis-eQTL analysis.")
@click.option('--maf_threshold', type=float, default=0.01, show_default=True, help="Minor allele frequency threshold.")
@click.option('--window', type=int, default=1000000, show_default=True, help="Window size (in base pairs) for cis-sQTL analysis.")
@click.option('--pval_threshold', type=float, default=1e-8, show_default=True, help="P-value threshold for trans-QTL analysis.")
def main(expression_bed, covariates_file, outfile, mode, nperm, maf_threshold, window, pval_threshold):
    """
    主函数，用于运行 QTL 分析。
    """
    # 加载表达数据和协变量
    phenotype_df, phenotype_pos_df = load_expression_data(expression_bed)
    covariates_df = load_covariates(covariates_file, phenotype_df.columns)

    # 加载基因型数据（使用硬编码路径）
    print(f"Loading PLINK data from: {PLINK_PREFIX_PATH}")
    pr = genotypeio.PlinkReader(PLINK_PREFIX_PATH, select_samples=phenotype_df.columns)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

    # 根据 mode 运行不同分析
    if mode == 'p':
        perform_cis_analysis(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, outfile, nperm, maf_threshold, window)
    elif mode == 'n':
        perform_nominal_mapping(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, outfile, maf_threshold, window)
    elif mode == 't':
        perform_trans_analysis(genotype_df, phenotype_df, covariates_df, outfile, pval_threshold, maf_threshold)

if __name__ == "__main__":
    main()