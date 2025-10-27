import click
import csv

@click.command()
@click.option('--cluster-file', required=True, help='Hap_gene_cluster.tsv')
@click.option('--so-expr-file', required=True, help='stage1_sum_expression_data_so_sums.tsv')
@click.option('--ss-expr-file', required=True, help='stage1_sum_expression_data_ss_sums.tsv')
@click.option('--output-file', default='hap_gene_expression_values.tsv', help='Output file name')
def extract(cluster_file, so_expr_file, ss_expr_file, output_file):
    """Extract gene-wise expression values from cluster-based expression matrices."""

    def load_expr(file_path):
        expr = {}
        with open(file_path, newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            cul_cols = [col for col in reader.fieldnames if col.startswith("Cul")]
            for row in reader:
                cid = f"cluster{row['Cluster']}"
                expr[cid] = {col: row[col] for col in cul_cols}
        return expr, cul_cols

    so_expr, cul_columns = load_expr(so_expr_file)
    ss_expr, _ = load_expr(ss_expr_file)

    records = []
    with open(cluster_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            cluster_id = f"cluster{row['Cluster']}"
            # 处理 so.Hap_genes
            so_genes = row.get("so.Hap_genes", "")
            if so_genes:
                for gene in so_genes.replace('，', ',').split(','):
                    gene = gene.strip()
                    if gene and cluster_id in so_expr:
                        expr = so_expr[cluster_id]
                        records.append([gene] + [expr.get(col, "0") for col in cul_columns])
            # 处理 ss.Hap_genes
            ss_genes = row.get("ss.Hap_genes", "")
            if ss_genes:
                for gene in ss_genes.replace('，', ',').split(','):
                    gene = gene.strip()
                    if gene and cluster_id in ss_expr:
                        expr = ss_expr[cluster_id]
                        records.append([gene] + [expr.get(col, "0") for col in cul_columns])

    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Gene"] + cul_columns)
        writer.writerows(records)

    click.echo(f"[✓] 输出完成：{output_file}")

if __name__ == '__main__':
    extract()