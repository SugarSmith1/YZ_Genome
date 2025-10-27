import pandas as pd
import click
import numpy as np
import re

@click.command()
@click.option('--cluster-file', required=True, help='Input TSV file with cluster gene lists')
@click.option('--expression-file', required=True, help='Input TSV file with expression data')
@click.option('--output-prefix', default='cluster', help='Prefix for output files')
@click.option('--chunk-size', default=1000, help='Number of clusters to process at once')
def process_clusters(cluster_file, expression_file, output_prefix, chunk_size):
    """Process clusters with prefix removal and expression-based deduplication."""
    
    # 1. 加载表达数据
    click.echo("Loading expression data...")
    expr_df = pd.read_csv(expression_file, sep='\t')
    
    if 'target_id' not in expr_df.columns:
        click.echo("ERROR: Expression data must contain 'target_id' column")
        return
    
    expr_df.set_index('target_id', inplace=True)
    expr_matrix = expr_df.values.astype(np.float32)
    gene_index = {str(gene): idx for idx, gene in enumerate(expr_df.index)}
    varieties = expr_df.columns.tolist()
    
    # 2. 初始化输出文件
    writers = {
        'so': open(f'{output_prefix}_so_sums.tsv', 'w'),
        'ss': open(f'{output_prefix}_ss_sums.tsv', 'w'),
        'combined': open(f'{output_prefix}_combined_sums.tsv', 'w')
    }
    
    # 写入表头
    headers = ['Cluster', 'so.Hap_genes', 'ss.Hap_genes'] + varieties
    for writer in writers.values():
        writer.write('\t'.join(headers) + '\n')
    
    # 3. 改进的ID标准化函数
    def normalize_gene_id(raw_id):
        """处理形如 so.YZ081609401660 的ID格式"""
        # 匹配前缀格式：so. 或 ss. 开头
        match = re.match(r'^(so|ss)\.(.*)$', raw_id)
        if match:
            return match.group(2)  # 返回前缀后的部分
        return raw_id  # 如果没有前缀则原样返回
    
    # 4. 基于表达量的去重函数
    def get_unique_vectors(genes):
        seen = set()
        unique_vectors = []
        for raw_gene in genes:
            # 进行ID标准化
            norm_gene = normalize_gene_id(raw_gene)
            
            # 匹配表达矩阵
            if norm_gene in gene_index:
                vec = expr_matrix[gene_index[norm_gene]]
                vec_rounded = tuple(np.round(vec, 4))
                if vec_rounded not in seen:
                    seen.add(vec_rounded)
                    unique_vectors.append(vec)
        return unique_vectors
    
    # 5. 处理每个cluster
    click.echo("Processing clusters with prefix removal...")
    for chunk in pd.read_csv(cluster_file, sep='\t', chunksize=chunk_size):
        for _, row in chunk.iterrows():
            cluster_id = row.name
            
            # 获取原始基因列表
            so_genes_raw = row['so.YZ081609_genes'] if pd.notna(row['so.YZ081609_genes']) else 'NA'
            ss_genes_raw = row['ss.YZ081609_genes'] if pd.notna(row['ss.YZ081609_genes']) else 'NA'
            
            # 分割基因列表
            so_genes = so_genes_raw.split('，') if so_genes_raw != 'NA' else []
            ss_genes = ss_genes_raw.split('，') if ss_genes_raw != 'NA' else []
            
            # 计算各组总和
            base_data = [f'cluster{cluster_id}', so_genes_raw, ss_genes_raw]
            
            # SO总和计算
            so_unique = get_unique_vectors(so_genes)
            so_sums = np.sum(so_unique, axis=0) if so_unique else np.zeros(len(varieties))
            writers['so'].write('\t'.join(map(str, base_data + [f"{x:.4f}" for x in so_sums])) + '\n')
            
            # SS总和计算
            ss_unique = get_unique_vectors(ss_genes)
            ss_sums = np.sum(ss_unique, axis=0) if ss_unique else np.zeros(len(varieties))
            writers['ss'].write('\t'.join(map(str, base_data + [f"{x:.4f}" for x in ss_sums])) + '\n')
            
            # 合并总和计算
            combined_genes = so_genes + ss_genes
            combined_unique = get_unique_vectors(combined_genes)
            combined_sums = np.sum(combined_unique, axis=0) if combined_unique else np.zeros(len(varieties))
            writers['combined'].write('\t'.join(map(str, base_data + [f"{x:.4f}" for x in combined_sums])) + '\n')
    
    # 关闭文件
    for writer in writers.values():
        writer.close()
    
    click.echo(f"Processing complete! Results saved to {output_prefix}_*_sums.tsv")

if __name__ == '__main__':
    process_clusters()
