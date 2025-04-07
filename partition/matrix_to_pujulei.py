import time
import numpy as np
from sklearn.cluster import SpectralClustering
from scipy.sparse.csgraph import laplacian

from Utils.read_qasm import list_str_to_int, converter_circ_from_qasm


def spectral_min_cut(graph, k):
    """
    使用谱聚类算法实现k个分区的最小割
    :param graph: 二维列表表示的图，graph[i][j]表示节点i到节点j的权重（边的频次）
    :param k: 想要分成的区数
    :return: 分区结果
    """

    start_time = time.time()

    # 将输入图转换为 NumPy 矩阵
    graph = np.array(graph)

    # 计算图的拉普拉斯矩阵
    L = laplacian(graph, normed=True)

    # 使用谱聚类算法
    clustering = SpectralClustering(n_clusters=k, affinity='precomputed', random_state=42)

    # 计算聚类结果
    labels = clustering.fit_predict(graph)

    print(labels)

    # 计算最小割的权重
    cut_weight = 0.0
    N = len(labels)
    for i in range(N):
        for j in range(i + 1, N):
            if labels[i] != labels[j]:
                cut_weight += graph[i][j]

    # print("最小割权重：", cut_weight)


    # end_time = time.time()
    # execution_time = end_time - start_time
    # print(execution_time)

    return labels


def gate_list_to_graph(gate_list):
    # 找到最大节点编号，决定邻接矩阵的大小
    max_node = max(max(pair) for pair in gate_list)

    # 初始化邻接矩阵，大小为 (max_node+1) x (max_node+1)
    graph = np.zeros((max_node + 1, max_node + 1), dtype=int)

    # 遍历 gate_list，累加边的频次
    for i, j in gate_list:
        graph[i][j] += 1

    return graph


if __name__ == '__main__':
    # # 示例图，二维列表表示，graph[i][j]表示节点i到节点j的权重（频次）
    # graph = [
    #     [0, 10, 1, 0],
    #     [10, 0, 2, 8],
    #     [1, 2, 0, 3],
    #     [0, 8, 3, 0]
    # ]

    # gate_list = [[0, 1], [2, 1], [1, 2], [3, 2], [2, 3], [4, 3], [3, 4], [4, 3], [4, 2], [4, 1], [4, 0], [0, 4], [0, 4],
    #              [1, 4], [1, 4], [2, 4], [2, 4], [1, 4], [1, 2]]
    # file_list = ['100', '200', '500', '1000', '2000', '5000', '10000']
    # file_list = ['200','300', '400']

    # file_list = ['2000', '3000', '4000']

    # file_list = ['500', '600', '700', '800', '900']
    # file_list = ['2000', '3000', '4000','5000']
    # # file_list = ['60', '70', '80', '90', '100']
    #
    # # file_list = ['41', '42', '43', '44', '45', '46', '47', '48', '49']
    # for file in file_list:
    #     input_filename = f'./random_qasm/{file}_10000.qasm'
    #
    #     gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    #     # 转换为邻接矩阵
    #     graph = gate_list_to_graph(gate_list)
    #     k = 50  # 想要分成的区数
    #
    #     result = spectral_min_cut(graph, k)

        # print("分区结果:", result)

    input_filename = f'../qasm/qasm_czl/mini_alu_305.qasm'

    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    graph = gate_list_to_graph(gate_list)

    spectral_min_cut(graph, 3)

