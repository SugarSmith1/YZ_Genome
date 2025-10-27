import click
import pandas as pd
import numpy as np
from scipy import stats
import logging

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

@click.command()
@click.option('--input', '-i', required=True, help='Input trans-eQTL file (TSV format)')
@click.option('--output', '-o', required=True, help='Output file path')
@click.option('--qval-threshold', default=0.05, show_default=True, help='FDR q-value threshold')
@click.option('--pval-threshold', default=1e-8, show_default=True, help='Raw p-value threshold (additional filter)')
@click.option('--effect-size-threshold', default=0.1, show_default=True, help='Minimum absolute effect size |slope|')
@click.option('--min-maf', default=0.05, show_default=True, help='Minimum minor allele frequency')
@click.option('--gene-count-threshold', default=3, show_default=True, help='Minimum number of genes per SNP (for hotspot detection)')
@click.option('--snps-per-gene', default=5, show_default=True, help='Maximum SNPs per gene to keep (top by p-value)')
@click.option('--filter-hotspots', is_flag=True, help='Filter out trans-eQTL hotspots')
@click.option('--filter-cis-acting', is_flag=True, help='Filter out cis-acting trans-eQTLs')
def filter_trans_eqtls(input, output, qval_threshold, pval_threshold, effect_size_threshold, 
                      min_maf, gene_count_threshold, snps_per_gene, filter_hotspots, filter_cis_acting):
    """
    过滤 trans-eQTL 结果，减少假阳性并提高结果质量
    """
    
    logger.info(f"Loading trans-eQTL data from: {input}")
    
    try:
        # 读取数据
        df = pd.read_csv(input, sep='\t')
        
        # 检查必要的列
        required_cols = ['phenotype_id', 'variant_id', 'pval', 'qval', 'slope']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logger.error(f"Missing required columns: {missing_cols}")
            return
        
        logger.info(f"Initial trans-eQTL count: {len(df)}")
        
        # 第一步：基本过滤
        filtered_df = df.copy()
        
        # FDR 过滤
        initial_count = len(filtered_df)
        filtered_df = filtered_df[filtered_df['qval'] < qval_threshold]
        logger.info(f"After FDR filter (qval < {qval_threshold}): {len(filtered_df)} / {initial_count}")
        
        # p-value 过滤
        initial_count = len(filtered_df)
        filtered_df = filtered_df[filtered_df['pval'] < pval_threshold]
        logger.info(f"After p-value filter (pval < {pval_threshold}): {len(filtered_df)} / {initial_count}")
        
        # 效应大小过滤
        initial_count = len(filtered_df)
        filtered_df = filtered_df[filtered_df['slope'].abs() > effect_size_threshold]
        logger.info(f"After effect size filter (|slope| > {effect_size_threshold}): {len(filtered_df)} / {initial_count}")
        
        # 第二步：MAF 过滤（如果有MAF列）
        if 'maf' in df.columns:
            initial_count = len(filtered_df)
            filtered_df = filtered_df[filtered_df['maf'] >= min_maf]
            logger.info(f"After MAF filter (MAF >= {min_maf}): {len(filtered_df)} / {initial_count}")
        
        # 第三步：限制每个基因的SNP数量
        initial_count = len(filtered_df)
        filtered_df = filtered_df.sort_values(['phenotype_id', 'pval'])
        filtered_df = filtered_df.groupby('phenotype_id').head(snps_per_gene)
        logger.info(f"After limiting to {snps_per_gene} SNPs per gene: {len(filtered_df)} / {initial_count}")
        
        # 第四步：过滤 trans-eQTL 热点（可选）
        if filter_hotspots:
            initial_count = len(filtered_df)
            # 计算每个SNP影响的基因数量
            snp_gene_counts = filtered_df.groupby('variant_id')['phenotype_id'].nunique()
            # 过滤掉影响过多基因的SNP（可能是技术假象）
            hotspot_snps = snp_gene_counts[snp_gene_counts >= gene_count_threshold].index
            filtered_df = filtered_df[~filtered_df['variant_id'].isin(hotspot_snps)]
            logger.info(f"After hotspot filter (SNPs affecting >={gene_count_threshold} genes): {len(filtered_df)} / {initial_count}")
        
        # 第五步：过滤 cis-acting trans-eQTL（可选）
        if filter_cis_acting and 'phenotype_chr' in df.columns and 'variant_chr' in df.columns:
            initial_count = len(filtered_df)
            # 识别染色体相同的关联（可能是cis泄漏）
            filtered_df = filtered_df[filtered_df['phenotype_chr'] != filtered_df['variant_chr']]
            logger.info(f"After removing cis-acting trans-eQTLs: {len(filtered_df)} / {initial_count}")
        
        # 第六步：计算统计量并排序
        if not filtered_df.empty:
            # 添加效应大小绝对值列
            filtered_df['abs_slope'] = filtered_df['slope'].abs()
            
            # 按p值排序
            filtered_df = filtered_df.sort_values(['pval', 'abs_slope'], ascending=[True, False])
            
            # 计算基本统计
            logger.info(f"Final trans-eQTL count: {len(filtered_df)}")
            logger.info(f"Unique genes: {filtered_df['phenotype_id'].nunique()}")
            logger.info(f"Unique SNPs: {filtered_df['variant_id'].nunique()}")
            logger.info(f"Median |effect size|: {filtered_df['abs_slope'].median():.4f}")
            logger.info(f"Min p-value: {filtered_df['pval'].min():.2e}")
            
            # 保存结果
            filtered_df.to_csv(output, sep='\t', index=False)
            logger.info(f"Filtered trans-eQTLs saved to: {output}")
            
            # 保存简要统计
            stats_output = output.replace('.txt', '_stats.txt')
            with open(stats_output, 'w') as f:
                f.write(f"Trans-eQTL Filtering Statistics\n")
                f.write(f"==============================\n")
                f.write(f"Input file: {input}\n")
                f.write(f"Output file: {output}\n")
                f.write(f"Final associations: {len(filtered_df)}\n")
                f.write(f"Unique genes: {filtered_df['phenotype_id'].nunique()}\n")
                f.write(f"Unique SNPs: {filtered_df['variant_id'].nunique()}\n")
                f.write(f"Median |effect size|: {filtered_df['abs_slope'].median():.4f}\n")
                f.write(f"Min p-value: {filtered_df['pval'].min():.2e}\n")
                f.write(f"Max p-value: {filtered_df['pval'].max():.2e}\n")
            
        else:
            logger.warning("No trans-eQTLs passed filtering criteria")
            # 创建空文件保持一致性
            pd.DataFrame(columns=df.columns).to_csv(output, sep='\t', index=False)
            
    except Exception as e:
        logger.error(f"Error processing trans-eQTL file: {e}")
        raise

if __name__ == '__main__':
    filter_trans_eqtls()