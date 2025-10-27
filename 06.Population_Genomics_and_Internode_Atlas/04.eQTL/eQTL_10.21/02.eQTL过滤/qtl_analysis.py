# -*- coding: utf-8 -*-
'''
Correct QTL analysis script - cis-eQTL outputs lead SNP per gene, trans-eQTL outputs all significant pairs
'''

import os
import torch
import click
import pandas as pd
import tensorqtl
from tensorqtl import genotypeio, cis, trans
from statsmodels.stats.multitest import multipletests

# 设置 CUDA
os.environ['CUDA_VISIBLE_DEVICES'] = "0"
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f'Using device: {device}')

class Config:
    PLINK_PREFIX_PATH = "/data0/agis_xiazhongqiang/Project/04.eQTL/06.YZ.qtl_mapping/01.data/07.pre.all.data/GWAS"
    CIS_NPERM = 10000
    CIS_WINDOW = 1000000
    MAF_THRESHOLD = 0.05
    FDR_THRESHOLD = 0.05
    TRANS_PVAL_THRESHOLD = 1e-5
    SEED = 2022

def load_data(expression_bed, covariates_file):
    """加载所有必需数据"""
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    print(f"Loaded expression data: {phenotype_df.shape[1]} samples, {phenotype_df.shape[0]} genes")
    
    # 加载协变量
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0)
    common_samples = list(set(phenotype_df.columns) & set(covariates_df.index))
    covariates_df = covariates_df.loc[common_samples]
    phenotype_df = phenotype_df[common_samples]
    print(f"Loaded covariates: {covariates_df.shape[1]} covariates for {len(common_samples)} samples")
    
    # 加载基因型
    print(f"Loading genotype data from: {Config.PLINK_PREFIX_PATH}")
    pr = genotypeio.PlinkReader(Config.PLINK_PREFIX_PATH, select_samples=phenotype_df.columns)
    genotype_df = pr.load_genotypes()
    variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]
    print(f"Loaded genotype data: {genotype_df.shape[1]} samples, {genotype_df.shape[0]} variants")
    
    return phenotype_df, phenotype_pos_df, covariates_df, genotype_df, variant_df

def run_cis_eqtl(genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df, output_prefix):
    """运行cis-eQTL分析，输出每个基因的lead SNP"""
    print("Running cis-eQTL analysis...")
    
    try:
        # cis映射 - 默认返回每个基因最显著的SNP
        cis_df = cis.map_cis(
            genotype_df, variant_df, phenotype_df, phenotype_pos_df, covariates_df,
            nperm=Config.CIS_NPERM, maf_threshold=Config.MAF_THRESHOLD, 
            window=Config.CIS_WINDOW, seed=Config.SEED
        )
        
        # FDR校正
        cis_df = cis.calculate_qvalues(cis_df, fdr=Config.FDR_THRESHOLD)
        
        # 过滤显著结果 - 每个基因已经是lead SNP
        significant_cis = cis_df[cis_df['qval'] < Config.FDR_THRESHOLD].copy()
        
        # 添加效应方向
        significant_cis['effect_direction'] = significant_cis['slope'].apply(
            lambda x: 'positive' if x > 0 else 'negative'
        )
        
        if not significant_cis.empty:
            output_file = f"{output_prefix}_cis_lead_snps.txt"
            significant_cis.to_csv(output_file, sep="\t", index=True)
            print(f"Cis-eQTL: {len(significant_cis)} significant genes -> {output_file}")
            return significant_cis
        else:
            print("Cis-eQTL: No significant associations")
            # 创建空文件保持一致性
            pd.DataFrame().to_csv(f"{output_prefix}_cis_lead_snps.txt", sep="\t")
            return pd.DataFrame()
            
    except Exception as e:
        print(f"Error in cis-eQTL analysis: {e}")
        # 创建空文件
        pd.DataFrame().to_csv(f"{output_prefix}_cis_lead_snps.txt", sep="\t")
        return pd.DataFrame()

def run_trans_eqtl(genotype_df, phenotype_df, covariates_df, output_prefix):
    """运行trans-eQTL分析，输出所有显著SNP-基因对"""
    print("Running trans-eQTL analysis...")
    
    try:
        # trans映射 - 返回所有达到阈值的SNP-基因对
        trans_df = trans.map_trans(
            genotype_df, phenotype_df, covariates_df,
            return_sparse=True, pval_threshold=Config.TRANS_PVAL_THRESHOLD,
            maf_threshold=Config.MAF_THRESHOLD, batch_size=10000
        )
        
        if trans_df.empty:
            print("Trans-eQTL: No associations found in initial screening")
            # 创建空文件
            pd.DataFrame().to_csv(f"{output_prefix}_trans_all_significant.txt", sep="\t")
            return pd.DataFrame()
        
        print(f"Trans-eQTL initial screening: {len(trans_df)} associations found")
        
        # FDR校正
        _, qvals, _, _ = multipletests(trans_df['pval'].values, method='fdr_bh')
        trans_df['qval'] = qvals
        
        # 过滤显著结果 - 保留所有显著对
        significant_trans = trans_df[trans_df['qval'] < Config.FDR_THRESHOLD].copy()
        
        if not significant_trans.empty:
            # 添加效应方向
            significant_trans['effect_direction'] = significant_trans['slope'].apply(
                lambda x: 'positive' if x > 0 else 'negative'
            )
            
            # 按p值排序
            significant_trans = significant_trans.sort_values('pval')
            
            output_file = f"{output_prefix}_trans_all_significant.txt"
            significant_trans.to_csv(output_file, sep="\t", index=False)
            print(f"Trans-eQTL: {len(significant_trans)} significant pairs -> {output_file}")
            return significant_trans
        else:
            print("Trans-eQTL: No significant associations after FDR correction")
            # 创建空文件
            pd.DataFrame().to_csv(f"{output_prefix}_trans_all_significant.txt", sep="\t")
            return pd.DataFrame()
            
    except Exception as e:
        print(f"Error in trans-eQTL analysis: {e}")
        # 创建空文件
        pd.DataFrame().to_csv(f"{output_prefix}_trans_all_significant.txt", sep="\t")
        return pd.DataFrame()

@click.command()
@click.option('--expression_bed', required=True, help="Expression BED file path")
@click.option('--covariates_file', required=True, help="Covariates file path")
@click.option('--outfile', required=True, help="Output file prefix")
@click.option('--mode', type=click.Choice(['p', 't', 'both']), required=True, 
              help="p=cis-eQTL, t=trans-eQTL, both=both analyses")
def main(expression_bed, covariates_file, outfile, mode):
    """
    QTL分析脚本:
    - cis-eQTL: 每个基因输出一个lead SNP
    - trans-eQTL: 输出所有显著SNP-基因对
    """
    print(f"Starting QTL analysis: mode={mode}")
    print(f"Expression file: {expression_bed}")
    print(f"Covariates file: {covariates_file}")
    print(f"Output prefix: {outfile}")
    
    # 加载数据
    phenotype_df, phenotype_pos_df, covariates_df, genotype_df, variant_df = load_data(
        expression_bed, covariates_file
    )
    
    # 运行指定分析
    if mode in ['p', 'both']:
        cis_results = run_cis_eqtl(
            genotype_df, variant_df, phenotype_df, phenotype_pos_df, 
            covariates_df, outfile
        )
    
    if mode in ['t', 'both']:
        trans_results = run_trans_eqtl(
            genotype_df, phenotype_df, covariates_df, outfile
        )
    
    print("QTL analysis completed successfully!")

if __name__ == "__main__":
    main()