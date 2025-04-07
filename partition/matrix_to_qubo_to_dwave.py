import numpy as np
import time
from Utils.read_qasm import list_str_to_int, converter_circ_from_qasm
import pandas as pd
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler, EmbeddingComposite, LeapHybridCQMSampler
import hybrid


def generate_upper_triangular_matrix(cnot_list, num_qubits):
    # 初始化一个n×n的零矩阵
    matrix = np.zeros((num_qubits, num_qubits), dtype=int)

    # 遍历CNOT门列表，更新上三角矩阵
    for control, target in cnot_list:
        if control < target:  # 仅更新上三角部分
            matrix[control, target] += 1

    return matrix


def generate_matrix(w, fai=0.5, lmd1=6, lmd2=5, rou=0):
    two_lmd = 2 * lmd1 + 2 * lmd2
    diagonal = lmd1 - 2 + lmd2 - 2 * (3 - rou)

    # 确定w的大小
    n = len(w)

    # 初始化一个nxn的matrix
    matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i, n):
            if i == j:
                matrix[i][j] = diagonal
            else:
                matrix[i][j] = fai - (1 + fai) * w[i][j] + two_lmd

    return -matrix




# def cut_circuit_by_qubo(matrix):
#     start_time = time.time()
#
#     matrix_n = kw.cim.normalizer(matrix, normalization=0.5)
#     output = kw.cim.simulator_core(
#         matrix_n,
#         c=0,
#         pump=0.7,
#         noise=0.01,
#         laps=1000,
#         dt=0.1)
#
#     h = kw.sampler.hamiltonian(matrix, output)
#
#     c_set = kw.sampler.binarizer(output)
#     opt = kw.sampler.optimal_sampler(matrix, c_set, 0)
#     best = opt[0][0]
#     # print(best)
#
#     end_time = time.time()
#     execution_time = end_time - start_time
#     print(execution_time)


# def generate_qubo_matrix_from_graph(graph, num_partitions):
#     num_nodes = len(graph)
#
#     # 初始化 QUBO 矩阵，大小为 (num_nodes * num_partitions, num_nodes * num_partitions)
#     qubo_matrix = np.zeros((num_nodes * num_partitions, num_nodes * num_partitions))
#
#     # 构建约束，确保每个节点只分配到一个分区
#     for i in range(num_nodes):
#         for k in range(num_partitions):
#             idx_ik = i * num_partitions + k
#             qubo_matrix[idx_ik, idx_ik] += 1  # 确保每个节点属于某个分区
#             for kp in range(num_partitions):
#                 if k != kp:
#                     idx_ikp = i * num_partitions + kp
#                     qubo_matrix[idx_ik, idx_ikp] += 2  # 惩罚节点被分配到多个分区
#
#     # 添加边权重作为目标函数的一部分
#     for i in range(num_nodes):
#         for j in range(i + 1, num_nodes):
#             weight = graph[i][j]
#             if weight != 0:  # 只考虑有边的情况
#                 for k in range(num_partitions):
#                     idx_ik = i * num_partitions + k
#                     idx_jk = j * num_partitions + k
#                     qubo_matrix[idx_ik, idx_jk] -= weight  # 如果 i 和 j 在同一分区，减少边的代价
#
#     return qubo_matrix


def generate_qubo_matrix_from_graph(graph, num_partitions):
    num_nodes = len(graph)

    # 初始化 QUBO 矩阵，大小为 (num_nodes * num_partitions, num_nodes * num_partitions)
    qubo_matrix = np.zeros((num_nodes * num_partitions, num_nodes * num_partitions))

    # 构建约束，确保每个节点只分配到一个分区
    for i in range(num_nodes):
        for k in range(num_partitions):
            idx_ik = i * num_partitions + k
            qubo_matrix[idx_ik, idx_ik] += 1  # 确保每个节点属于某个分区
            for kp in range(k + 1, num_partitions):  # 确保只填充上三角部分
                idx_ikp = i * num_partitions + kp
                qubo_matrix[idx_ik, idx_ikp] += 2  # 惩罚节点被分配到多个分区

    # 添加边权重作为目标函数的一部分
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):  # 确保只填充上三角部分
            weight = graph[i][j]
            if weight != 0:  # 只考虑有边的情况
                for k in range(num_partitions):
                    idx_ik = i * num_partitions + k
                    idx_jk = j * num_partitions + k
                    qubo_matrix[idx_ik, idx_jk] -= weight  # 如果 i 和 j 在同一分区，减少边的代价

    # 调整惩罚系数，避免所有节点分配到一个分区
    for i in range(num_nodes):
        for k in range(1, num_partitions):
            idx_ik = i * num_partitions + k
            qubo_matrix[idx_ik, idx_ik] += 1.9  # 降低多分区的惩罚系数，以鼓励分配到不同分区

    # 添加惩罚，确保每个分区至少有两个节点
    for k in range(num_partitions):
        for i in range(num_nodes):
            for j in range(i + 1, num_nodes):  # 检查每一对节点是否属于同一分区
                idx_ik = i * num_partitions + k
                idx_jk = j * num_partitions + k
                if graph[i][j] != 0:  # 如果节点 i 和 j 之间有边
                    qubo_matrix[idx_ik, idx_jk] += 50  # 增加惩罚，确保每个分区至少有两个节点

    # 进一步加大对于割切边的惩罚，避免更多的跨分区连接
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            weight = graph[i][j]
            if weight != 0:
                for k in range(num_partitions):
                    idx_ik = i * num_partitions + k
                    idx_jk = j * num_partitions + k
                    # 惩罚不同分区的连接
                    qubo_matrix[idx_ik, idx_jk] += weight  # 增加跨分区的边权重惩罚

    return qubo_matrix




if __name__ == '__main__':
        # 示例CNOT门列表

        input_filename = f'../qasm/4gt5_76.qasm'
        gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
        num_qubits = max([max(row) for row in gate_list]) + 1  # 量子比特数量

        start_time = time.time()

        # 生成上三角矩阵 表示图
        upper_triangular_matrix = generate_upper_triangular_matrix(gate_list, num_qubits)

        # print(upper_triangular_matrix)

        k = 2 # 分区数

        qubo_matrix = generate_qubo_matrix_from_graph(upper_triangular_matrix, k)
        # qubo_matrix= np.random.randint(0, 50, (10, 10))
        # qubo_matrix = generate_matrix2()
        print(qubo_matrix)

        # print(len(qubo_matrix))

        # # 调用函数将 NumPy 矩阵转换为 QUBO 字典
        # qubo = numpy_to_qubo(qubo_matrix)
        # print("QUBO 字典: ", qubo)
        #
        # # 创建一个二进制二次模型 (Binary Quadratic Model)
        # bqm = BinaryQuadraticModel.from_qubo(qubo)
        #
        # # 打印 Binary Quadratic Model
        # print("Binary Quadratic Model: ", bqm)



        # dwave真机
        bqm = BinaryQuadraticModel.from_numpy_matrix(qubo_matrix)
        # 打印 Binary Quadratic Model
        print("Binary Quadratic Model: ", bqm)

        # # 使用 D-Wave 系统的量子退火机
        # sampler = EmbeddingComposite(DWaveSampler())
        # # 提交任务并获取结果
        # sampleset = sampler.sample(bqm, num_reads=100, label='file')
        # # 打印结果
        # print(sampleset)
        #
        # end_time1 = time.time()
        # execution_time1 = end_time1 - start_time


        # 使用 D-Wave 系统的hybird混合量子退火机
        # Define the workflow
        iteration = hybrid.RacingBranches(
            hybrid.InterruptableTabuSampler(),
            hybrid.EnergyImpactDecomposer(size=2)
            | hybrid.QPUSubproblemAutoEmbeddingSampler()
            | hybrid.SplatComposer()
        ) | hybrid.ArgMin()
        workflow = hybrid.LoopUntilNoImprovement(iteration, convergence=3)

        # Solve the problem
        init_state = hybrid.State.from_problem(bqm)
        final_state = workflow.run(init_state).result()


        # Print results
        print("Solution: sample={.samples}".format(final_state))

        # end_time = time.time()
        # execution_time = end_time - start_time
        # print( execution_time1)
        # print('总时间',execution_time)
        # print('与云平台通信时间', execution_time-execution_time1)


        # # 将 QUBO 矩阵保存为 CSV 文件
        # file_path_csv = f'./CSV/{file}_1000.csv'
        # df = pd.DataFrame(qubo_matrix)
        # df.to_csv(file_path_csv, index=False,header=False)
        #
        #
        # # matrix = generate_matrix(upper_triangular_matrix)
        # # print(matrix)
        # #
        # cut_circuit_by_qubo(qubo_matrix)