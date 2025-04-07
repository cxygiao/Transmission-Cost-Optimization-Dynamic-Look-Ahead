import numpy as np
import time
import kaiwu as kw
from Utils.read_qasm import list_str_to_int, converter_circ_from_qasm


def generate_upper_triangular_matrix(cnot_list, num_qubits):
    # 初始化一个n×n的零矩阵
    matrix = np.zeros((num_qubits, num_qubits), dtype=int)

    # 遍历CNOT门列表，更新上三角矩阵
    for control, target in cnot_list:
        if control < target:  # 仅更新上三角部分
            matrix[control, target] += 1

    return matrix


def generate_matrix(w, fai=0.5, lmd1=125, lmd2=86, rou=0):
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



def generate_matrix2():
    fai = 0.5
    lmd1 = 125
    lmd2 = 86
    rou = 0
    two_lmd = 2*lmd1+2*lmd2
    diagonal = lmd1-2 + lmd2-2*(3-rou)

    w=[[0,0,1,1,1,1],
       [0,0,0,1,2,0],
       [0,0,0,0,0,1],
       [0,0,0,0,0,0],
       [0,0,0,0,0,0],
       [0,0,0,0,0,0]
       ]

    matrix = -np.array([
        [diagonal, fai-(1+fai)*w[0][1]+two_lmd, fai-(1+fai)*w[0][2]+two_lmd, fai-(1+fai)*w[0][3]+two_lmd, fai-(1+fai)*w[0][4]+two_lmd, fai-(1+fai)*w[0][5]+two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd],
        [0, diagonal, fai-(1+fai)*w[1][2]+two_lmd, fai-(1+fai)*w[1][3]+two_lmd, fai-(1+fai)*w[1][4]+two_lmd, fai-(1+fai)*w[1][5]+two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd],
        [0, 0, diagonal, fai-(1+fai)*w[2][3]+two_lmd, fai-(1+fai)*w[2][4]+two_lmd, fai-(1+fai)*w[2][5]+two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd],
        [0, 0, 0, diagonal, fai-(1+fai)*w[3][4]+two_lmd, fai-(1+fai)*w[3][5]+two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd],
        [0, 0, 0, 0, diagonal, fai-(1+fai)*w[4][5]+two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd],
        [0, 0, 0, 0, 0, diagonal, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd, two_lmd],
        [0, 0, 0, 0, 0, 0, diagonal, fai-(1+fai)*w[0][1]+two_lmd, fai-(1+fai)*w[0][2]+two_lmd, fai-(1+fai)*w[0][3]+two_lmd, fai-(1+fai)*w[0][4]+two_lmd, fai-(1+fai)*w[0][5]+two_lmd],
        [0, 0, 0, 0, 0, 0, 0, diagonal, fai-(1+fai)*w[1][2]+two_lmd, fai-(1+fai)*w[1][3]+two_lmd, fai-(1+fai)*w[1][4]+two_lmd, fai-(1+fai)*w[1][5]+two_lmd],
        [0, 0, 0, 0, 0, 0, 0, 0, diagonal, fai-(1+fai)*w[2][3]+two_lmd, fai-(1+fai)*w[2][4]+two_lmd, fai-(1+fai)*w[2][5]+two_lmd],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, diagonal, fai-(1+fai)*w[3][4]+two_lmd, fai-(1+fai)*w[3][5]+two_lmd],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, diagonal, fai-(1+fai)*w[4][5]+two_lmd],
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, diagonal]
    ])
    return matrix

def generate_matrix_new(w):
    fai = 0.5
    lmd1 = 125
    lmd2 = 86
    rou = 0
    two_lmd = 2 * lmd1 + 2 * lmd2
    diagonal = lmd1 - 2 + lmd2 - 2 * (3 - rou)

    n = len(w)

    M = np.zeros((2 * n, 2 * n))

    for i in range(n):
        for j in range(n):
            if j < i:
                M[i, j] = 0
            elif j == i:
                M[i, j] = diagonal
            else:
                M[i, j] = fai - (1 + fai) * w[i][j] + two_lmd
        for j in range(n, 2 * n):
            M[i, j] = two_lmd

    for i in range(n, 2 * n):
        k = i - n
        for j in range(n):
            M[i, j] = 0
        for j in range(n, 2 * n):
            jprime = j - n
            if jprime < k:
                M[i, j] = 0
            elif jprime == k:
                M[i, j] = diagonal
            else:
                M[i, j] = fai - (1 + fai) * w[k][jprime] + two_lmd

    return -M



def cut_circuit_by_qubo(matrix):
    start_time = time.time()

    matrix_n = kw.cim.normalizer(matrix, normalization=0.5)
    output = kw.cim.simulator_core(
        matrix_n,
        c=0,
        pump=0.7,
        noise=0.01,
        laps=1000,
        dt=0.1)

    h = kw.sampler.hamiltonian(matrix, output)

    c_set = kw.sampler.binarizer(output)
    opt = kw.sampler.optimal_sampler(matrix, c_set, 0)
    best = opt[0][0]
    print(best)

    end_time = time.time()
    execution_time = end_time - start_time
    print(execution_time)


def generate_qubo_matrix_from_graph(graph, num_partitions):
    num_nodes = len(graph)

    qubo_matrix = np.zeros((num_nodes * num_partitions, num_nodes * num_partitions))

    # 构建约束，确保每个节点只分配到一个分区
    for i in range(num_nodes):
        for k in range(num_partitions):
            idx_ik = i * num_partitions + k
            qubo_matrix[idx_ik, idx_ik] += 1  # 确保每个节点属于某个分区
            for kp in range(num_partitions):
                if k != kp:
                    idx_ikp = i * num_partitions + kp
                    qubo_matrix[idx_ik, idx_ikp] += 2

    # 添加边权重作为目标函数的一部分
    for i in range(num_nodes):
        for j in range(i + 1, num_nodes):
            weight = graph[i][j]
            if weight != 0:  # 只考虑有边的情况
                for k in range(num_partitions):
                    idx_ik = i * num_partitions + k
                    idx_jk = j * num_partitions + k
                    qubo_matrix[idx_ik, idx_jk] -= weight

    return qubo_matrix


if __name__ == '__main__':
    # 示例CNOT门列表
    # gate_list = [[5, 43], [25, 43], [24, 35], [25, 35], [43, 35], [24, 35], [43, 35], [24, 35], [24, 43], [24, 43], [24, 35], [25, 35], [24, 35], [25, 35], [25, 24], [25, 24], [40, 26], [25, 26], [35, 26], [40, 26], [35, 26], [40, 26], [40, 35], [40, 35], [40, 26], [25, 26], [40, 26], [25, 26], [25, 40], [25, 40], [39, 1], [25, 1], [26, 1], [39, 1], [26, 1], [39, 1], [39, 26], [39, 26], [39, 1], [25, 1], [39, 1], [25, 1], [25, 39], [25, 39], [42, 1], [32, 1], [15, 1], [8, 1], [14, 1], [2, 1], [39, 17], [26, 17], [39, 17], [25, 17], [39, 17], [25, 17], [26, 17], [25, 39], [25, 39], [39, 17], [26, 17], [39, 17], [39, 26], [39, 26], [34, 17], [42, 17], [6, 17], [11, 17], [14, 17], [2, 17], [40, 41], [35, 41], [40, 41], [25, 41], [40, 41], [25, 41], [25, 40], [25, 40], [35, 41], [40, 41], [35, 41], [40, 41], [40, 35], [40, 35], [25, 9], [41, 9], [39, 9], [41, 9], [39, 9], [39, 41], [39, 41], [39, 9], [25, 9], [39, 9], [25, 9], [25, 39], [25, 39], [34, 9], [32, 9], [8, 9], [11, 9], [14, 9], [38, 9], [41, 44], [39, 44], [25, 44], [39, 44], [25, 44], [25, 39], [25, 39], [41, 44], [39, 44], [41, 44], [39, 44], [39, 41], [39, 41], [6, 44], [15, 44], [8, 44], [11, 44], [2, 44], [38, 44], [24, 18], [5, 18], [24, 18], [5, 18], [5, 24], [5, 24], [24, 18], [18, 23], [40, 23], [18, 23], [40, 23], [40, 18], [40, 18], [40, 23], [23, 36], [39, 36], [23, 36], [39, 36], [39, 23], [39, 23], [39, 36], [15, 36], [39, 26], [23, 26], [39, 26], [23, 26], [39, 26], [39, 23], [39, 23], [23, 26], [2, 26], [25, 27], [43, 27], [24, 27], [43, 27], [24, 27], [24, 43], [24, 43], [24, 27], [25, 27], [24, 27], [25, 27], [25, 24], [25, 24], [27, 28], [40, 28], [25, 28], [40, 28], [25, 28], [25, 40], [25, 40], [27, 28], [40, 28], [27, 28], [40, 28], [40, 27], [40, 27], [25, 29], [28, 29], [39, 29], [28, 29], [39, 29], [39, 28], [39, 28], [39, 29], [25, 29], [39, 29], [25, 29], [25, 39], [25, 39], [40, 30], [18, 30], [40, 30], [18, 30], [40, 30], [40, 18], [40, 18], [18, 30], [39, 31], [30, 31], [39, 31], [30, 31], [39, 31], [39, 30], [39, 30], [30, 31], [14, 31], [24, 4], [5, 4], [24, 4], [5, 4], [5, 24], [5, 24], [5, 4], [4, 16], [40, 16], [4, 16], [40, 16], [40, 4], [40, 4], [40, 16], [16, 0], [39, 0], [16, 0], [39, 0], [39, 16], [39, 16], [39, 0], [11, 0], [39, 12], [16, 12], [39, 12], [16, 12], [39, 12], [39, 16], [39, 16], [16, 12], [32, 12], [40, 19], [4, 19], [40, 19], [4, 19], [40, 19], [40, 4], [40, 4], [4, 19], [19, 3], [39, 3], [19, 3], [39, 3], [39, 19], [39, 19], [39, 3], [6, 3], [39, 7], [19, 7], [39, 7], [19, 7], [39, 7], [39, 19], [39, 19], [19, 7], [8, 7], [24, 10], [24, 10], [5, 10], [24, 10], [5, 10], [5, 24], [5, 24], [5, 10], [10, 33], [40, 33], [10, 33], [40, 33], [40, 10], [40, 10], [40, 33], [33, 29], [39, 29], [33, 29], [39, 29], [39, 33], [39, 33], [39, 29], [42, 29], [39, 30], [33, 30], [39, 30], [33, 30], [39, 30], [39, 33], [39, 33], [33, 30], [34, 30], [40, 28], [10, 28], [40, 28], [10, 28], [40, 28], [40, 10], [40, 10], [10, 28], [28, 22], [39, 22], [28, 22], [39, 22], [39, 28], [39, 28], [39, 22], [38, 22]]
    input_filename = './random_qasm/20_10000.qasm'
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))

    num_qubits = max([max(row) for row in gate_list]) + 1  # 量子比特数量

    # 生成上三角矩阵
    upper_triangular_matrix = generate_upper_triangular_matrix(gate_list, num_qubits)

    print(upper_triangular_matrix)

    qubo_matrix = generate_qubo_matrix_from_graph(upper_triangular_matrix, 2)
    print(qubo_matrix)

    # matrix = generate_matrix(upper_triangular_matrix)
    # print(matrix)
    #
    cut_circuit_by_qubo(qubo_matrix)