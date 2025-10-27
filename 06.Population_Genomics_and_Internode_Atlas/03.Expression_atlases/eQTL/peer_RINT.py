import numpy as np
import pandas as pd
import scipy.stats as ss
import click

def rank_INT(series, c=0.5, stochastic=False):
    """
    对 pandas Series 进行基于秩的逆正态变换。
    如果 stochastic 为 True，则对相同值的秩进行随机处理；否则，相同值共享相同的秩。NaN 值将被忽略。

    参数:
        series (pd.Series): 需要变换的 Series。
        c (float): Blom 常数，默认值为 0.5。
        stochastic (bool): 是否对相同值的秩进行随机处理，默认值为 False。

    返回:
        pd.Series: 变换后的 Series。
    """
    # 检查输入
    assert isinstance(series, pd.Series), "输入必须为 pandas Series"
    assert isinstance(c, float), "c 必须为浮点数"
    assert isinstance(stochastic, bool), "stochastic 必须为布尔值"

    # 设置随机种子
    np.random.seed(123)

    # 保存原始索引
    orig_idx = series.index

    # 去除 NaN 值
    series = series.loc[~pd.isnull(series)]

    # 计算秩
    if stochastic:
        # 随机打乱索引
        series = series.loc[np.random.permutation(series.index)]
        # 使用顺序法计算秩（相同值按位置区分）
        rank = ss.rankdata(series, method="ordinal")
    else:
        # 使用平均法计算秩（相同值共享平均秩）
        rank = ss.rankdata(series, method="average")

    # 将秩转换为 Series
    rank = pd.Series(rank, index=series.index)

    # 将秩转换为正态分布
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank))

    # 返回原始索引对应的结果
    return transformed[orig_idx]

def rank_to_normal(rank, c, n):
    """
    将秩转换为正态分布的分位数。

    参数:
        rank (float): 秩值。
        c (float): Blom 常数。
        n (int): 总样本数。

    返回:
        float: 正态分布的分位数。
    """
    # 计算标准分位数
    x = (rank - c) / (n - 2 * c + 1)
    return ss.norm.ppf(x)

@click.command()
@click.option("--peer_file", type=click.Path(exists=True), required=True, help="PEER 结果文件路径（CSV格式，行为样本，列为基因）。")
@click.option("--expre_file", type=click.Path(exists=True), required=True, help="expre结果文件路径（CSV格式，行为样本，列为基因）。")
@click.option("--output", type=click.Path(), required=True, help="输出文件路径（逆正态变换后的数据）。")
@click.option("--stochastic", is_flag=True, help="是否对相同值的秩进行随机处理。")
def main(peer_file, expre_file, output, stochastic):
    """
    对 PEER 结果文件进行逆正态变换，并保存结果。
    """
    # 读取 PEER 结果文件，并转置数据（行是基因，列是样本）
    peer_data = pd.read_csv(peer_file, header=None, index_col=0, sep="\t")
    peer_data = peer_data.T  # 转置，使得行是样本，列是基因

    # 读取表达数据文件
    SampleList = pd.read_csv(expre_file, header=0, index_col=0, sep="\t")

    # 对每列（基因）进行逆正态变换
    peer_data_transformed = peer_data.apply(lambda x: rank_INT(x, stochastic=stochastic), axis=0)

    # 修改索引和列，匹配表达数据的行名和列名
    peer_data_transformed.index = SampleList.columns  # 样本ID
    peer_data_transformed.columns = SampleList.index  # 基因名

    # 设置行和列名
    peer_data_transformed.index.name = 'IID'

    # 保存逆正态变换后的数据
    peer_data_transformed.to_csv(output, sep="\t", header=True, index=True)
    print(f"逆正态变换后的数据已保存至: {output}")

if __name__ == "__main__":
    main()