import click
import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_variant_id(variant_id):
    """解析变异ID，提取染色体和位置"""
    try:
        # 假设格式为 "chr1_100500_A_G" 或 "chr1:100500"
        if '_' in variant_id:
            chrom_part, pos_part = variant_id.split('_')[:2]
        elif ':' in variant_id:
            chrom_part, pos_part = variant_id.split(':')
        else:
            return None, None
        
        # 提取染色体
        chrom = chrom_part.replace('chr', '')
        # 提取位置
        pos = int(pos_part)
        return chrom, pos
    except:
        return None, None

def parse_gene_id(gene_id, gene_pos_dict):
    """从基因ID获取基因位置"""
    if gene_id in gene_pos_dict:
        return gene_pos_dict[gene_id]
    return None, None

@click.command()
@click.option('--input', '-i', required=True, help='Input trans-eQTL file')
@click.option('--output', '-o', required=True, help='Output file path')
@click.option('--gene-bed', '-g', required=True, help='Gene BED file with positions')
@click.option('--distance-threshold', default=5000000, show_default=True, 
              help='Minimum distance for same-chromosome trans-eQTL (bp)')
@click.option('--qval-threshold', default=0.05, show_default=True)
@click.option('--pval-threshold', default=1e-8, show_default=True)
@click.option('--effect-size-threshold', default=0.1, show_default=True)
@click.option('--snps-per-gene', default=5, show_default=True)
def filter_true_trans_eqtls(input, output, gene_bed, distance_threshold, 
                           qval_threshold, pval_threshold, effect_size_threshold, 
                           snps_per_gene):
    """
    过滤真正的 trans-eQTL：
    - 不同染色体，或
    - 同一染色体但距离 > 5Mb
    """
    
    logger.info("Loading gene position data...")
    # 加载基因位置信息
    gene_pos_df = pd.read_csv(gene_bed, sep='\t', header=None, 
                             names=['chrom', 'start', 'end', 'gene_id'])
    gene_pos_dict = {}
    for _, row in gene_pos_df.iterrows():
        gene_pos_dict[row['gene_id']] = (row['chrom'].replace('chr', ''), 
                                        (row['start'] + row['end']) // 2)  # 使用TSS
    
    logger.info(f"Loaded {len(gene_pos_dict)} gene positions")
    
    # 加载 trans-eQTL 结果
    logger.info(f"Loading trans-eQTL data from: {input}")
    df = pd.read_csv(input, sep='\t')
    
    if df.empty:
        logger.warning("Input file is empty")
        pd.DataFrame().to_csv(output, sep='\t', index=False)
        return
    
    logger.info(f"Initial trans-eQTL count: {len(df)}")
    
    # 第一步：基本统计过滤
    filtered_df = df.copy()
    
    # FDR 过滤
    filtered_df = filtered_df[filtered_df['qval'] < qval_threshold]
    logger.info(f"After FDR filter: {len(filtered_df)}")
    
    # p-value 过滤
    filtered_df = filtered_df[filtered_df['pval'] < pval_threshold]
    logger.info(f"After p-value filter: {len(filtered_df)}")
    
    # 效应大小过滤
    filtered_df = filtered_df[filtered_df['slope'].abs() > effect_size_threshold]
    logger.info(f"After effect size filter: {len(filtered_df)}")
    
    if filtered_df.empty:
        logger.warning("No associations passed basic filters")
        pd.DataFrame().to_csv(output, sep='\t', index=False)
        return
    
    # 第二步：识别真正的 trans-eQTL
    logger.info("Identifying true trans-eQTLs...")
    
    true_trans_eqtls = []
    
    for _, row in filtered_df.iterrows():
        gene_id = row['phenotype_id']
        variant_id = row['variant_id']
        
        # 获取基因位置
        gene_chrom, gene_pos = parse_gene_id(gene_id, gene_pos_dict)
        if gene_chrom is None:
            continue
        
        # 获取变异位置
        var_chrom, var_pos = parse_variant_id(variant_id)
        if var_chrom is None:
            continue
        
        # 判断是否为真正的 trans-eQTL
        is_true_trans = False
        trans_type = "unknown"
        
        if gene_chrom != var_chrom:
            # 不同染色体 - 肯定是 trans
            is_true_trans = True
            trans_type = "different_chrom"
        else:
            # 同一染色体，检查距离
            distance = abs(gene_pos - var_pos)
            if distance > distance_threshold:
                is_true_trans = True
                trans_type = f"same_chrom_{distance}bp"
        
        if is_true_trans:
            new_row = row.copy()
            new_row['gene_chrom'] = gene_chrom
            new_row['gene_pos'] = gene_pos
            new_row['var_chrom'] = var_chrom
            new_row['var_pos'] = var_pos
            new_row['distance'] = abs(gene_pos - var_pos) if gene_chrom == var_chrom else float('inf')
            new_row['trans_type'] = trans_type
            true_trans_eqtls.append(new_row)
    
    # 转换为DataFrame
    if true_trans_eqtls:
        true_trans_df = pd.DataFrame(true_trans_eqtls)
        logger.info(f"True trans-eQTLs identified: {len(true_trans_df)}")
        
        # 按染色体对类型分类
        chrom_types = true_trans_df['trans_type'].value_counts()
        logger.info("Trans-eQTL types:")
        for trans_type, count in chrom_types.items():
            logger.info(f"  {trans_type}: {count}")
        
        # 限制每个基因的SNP数量
        true_trans_df = true_trans_df.sort_values(['phenotype_id', 'pval'])
        true_trans_df = true_trans_df.groupby('phenotype_id').head(snps_per_gene)
        logger.info(f"After limiting to {snps_per_gene} SNPs per gene: {len(true_trans_df)}")
        
        # 排序并保存
        true_trans_df = true_trans_df.sort_values(['pval', 'distance'], ascending=[True, False])
        true_trans_df.to_csv(output, sep='\t', index=False)
        logger.info(f"True trans-eQTLs saved to: {output}")
        
        # 保存统计信息
        stats_file = output.replace('.txt', '_stats.txt')
        with open(stats_file, 'w') as f:
            f.write("True Trans-eQTL Statistics\n")
            f.write("==========================\n")
            f.write(f"Total true trans-eQTLs: {len(true_trans_df)}\n")
            f.write(f"Different chromosome: {len(true_trans_df[true_trans_df['trans_type'] == 'different_chrom'])}\n")
            f.write(f"Same chromosome >5Mb: {len(true_trans_df[true_trans_df['trans_type'].str.startswith('same_chrom')])}\n")
            f.write(f"Unique genes: {true_trans_df['phenotype_id'].nunique()}\n")
            f.write(f"Unique SNPs: {true_trans_df['variant_id'].nunique()}\n")
            f.write(f"Median distance (same chrom): {true_trans_df[true_trans_df['distance'] < float('inf')]['distance'].median():,.0f} bp\n")
        
    else:
        logger.warning("No true trans-eQTLs found")
        pd.DataFrame().to_csv(output, sep='\t', index=False)

if __name__ == '__main__':
    filter_true_trans_eqtls()