import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import click
import os

@click.command()
@click.option("--input1", type=click.Path(exists=True), required=True, help="第一个时期的表达数据文件路径（TSV格式，行为基因，列为样本）。")
@click.option("--input2", type=click.Path(exists=True), required=True, help="第二个时期的表达数据文件路径（TSV格式，行为基因，列为样本）。")
@click.option("--input3", type=click.Path(exists=True), required=True, help="第三个时期的表达数据文件路径（TSV格式，行为基因，列为样本）。")
@click.option("--input4", type=click.Path(exists=True), required=True, help="第四个时期的表达数据文件路径（TSV格式，行为基因，列为样本）。")
@click.option("--output_dir", type=click.Path(), required=True, help="输出目录路径，用于保存PCA结果文件和PCA图。")
def main(input1, input2, input3, input4, output_dir):
    """
    对四个时期的表达数据进行合并标准化，然后进行PCA分析，并绘制三维PCA图，每个时期用不同颜色区分。
    """
    # 定义颜色和标签
    colors = ['red', 'blue', 'green', 'purple']
    labels = ['Period1', 'Period2', 'Period3', 'Period4']
    inputs = [input1, input2, input3, input4]

    # 创建输出目录
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 读取所有文件并提取共同的基因集
    datasets = []
    for input_file in inputs:
        data = pd.read_csv(input_file, sep="\t", index_col=0)
        # 去掉有缺失值的基因行
        data = data.dropna()
        datasets.append(data)
        print(f"已读取文件: {input_file}, 样本数: {data.shape[1]}, 基因数: {data.shape[0]}")

    # 找到所有数据集中共同的基因
    common_genes = set(datasets[0].index)
    for data in datasets[1:]:
        common_genes.intersection_update(data.index)
    print(f"共同基因数: {len(common_genes)}")

    # 将共同基因转换为列表并排序
    common_genes = sorted(common_genes)

    # 合并所有时期的数据
    combined_data = pd.concat([data.loc[common_genes] for data in datasets], axis=1)
    print(f"合并后的数据, 样本数: {combined_data.shape[1]}, 基因数: {combined_data.shape[0]}")

    # 转置数据，使行为样本，列为基因
    combined_data = combined_data.T

    # 去除方差为0的基因
    combined_data = combined_data.loc[:, combined_data.std() > 0]
    print(f"去除方差为0的基因后, 基因数: {combined_data.shape[1]}")

    # 标准化数据（均值为0，方差为1）
    combined_data_standardized = (combined_data - combined_data.mean()) / combined_data.std()

    # 检查是否存在NaN
    if combined_data_standardized.isnull().any().any():
        print("标准化后存在NaN，将去除包含NaN的样本或基因。")
        combined_data_standardized = combined_data_standardized.dropna(axis=0)  # 去除包含NaN的样本
        combined_data_standardized = combined_data_standardized.dropna(axis=1)  # 去除包含NaN的基因

    # 执行PCA
    pca = PCA(n_components=3)
    pca_result = pca.fit_transform(combined_data_standardized)

    # 将PCA结果保存为DataFrame
    pca_df = pd.DataFrame(data=pca_result, columns=[f"PC{i+1}" for i in range(3)])
    pca_df.index = combined_data.index  # 保留样本名称

    # 为每个样本添加时期标签
    sample_periods = []
    for i, input_file in enumerate(inputs):
        data = pd.read_csv(input_file, sep="\t", index_col=0)
        sample_periods.extend([labels[i]] * data.shape[1])
    pca_df['Period'] = sample_periods

    # 保存PCA结果
    output_file = os.path.join(output_dir, "combined_pca_results.tsv")
    pca_df.to_csv(output_file, sep="\t")
    print(f"合并后的PCA结果已保存至: {output_file}")

    # 绘制三维PCA图
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # 为每个时期绘制散点图
    for i, label in enumerate(labels):
        period_data = pca_df[pca_df['Period'] == label]
        ax.scatter(period_data['PC1'], period_data['PC2'], period_data['PC3'], 
                   s=50, c=colors[i], label=label)

    # 设置坐标轴标签
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_zlabel("PC3")

    # 设置标题和图例
    ax.set_title("3D PCA Plot of Four Periods (Combined Standardization)")
    ax.legend()

    # 保存为PDF
    plot_file = os.path.join(output_dir, "combined_pca_plot.pdf")
    plt.savefig(plot_file, format="pdf")
    print(f"合并后的三维PCA图已保存至: {plot_file}")

if __name__ == "__main__":
    main()