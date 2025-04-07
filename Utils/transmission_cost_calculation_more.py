'''
功能：根据依赖图计算分布式线路的传输代价。
注意：此功能适用于多个分区的分布式线路的传输代价计算。
开始时间：2023.2.22
作者：陈新宇
版本: 2.0
完成时间：2023.2.26
'''

import copy
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, AncillaRegister, QuantumCircuit
import random

from Utils.read_qasm import get_data
from Utils.transmission_cost_calculation_2 import random_perturbations_transfer_qubit_list


Is_show_DAG = 0

# 随机扰动迭代次数
number_of_disturbance_iterations = 100

# 禁忌搜索参数设置
number_of_iterations = 50  # 迭代次数
number_of_times_not_updated = 5  # 最优解不更新次数 达到后直接返回禁忌表中最优解
length_of_taboo_table = 50  # 禁忌表长度
number_of_candidate_generation = 100  # 候选解数量


# 线路对应的依赖图
class Graph():
    def __init__(self, num_of_nodes, directed=True):
        self.num_of_nodes = num_of_nodes
        self.directed = directed
        self.edge_matrix = [[0 for _ in range(num_of_nodes)] for _ in range(num_of_nodes)]

    def add_edge(self, node1, node2, weight):
        self.edge_matrix[node1][node2] = weight
        if not self.directed:
            self.edge_matrix[node2][node1] = weight

    def print_edge_matrix(self):
        for i in range(self.num_of_nodes):
            print('G' + str(i), end=' ')
            print(self.edge_matrix[i])

# 绘制线路
def draw_new_circuit(gate_list):
    # 量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    qubits = QuantumRegister(circuit_qubit, 'q')
    # 生成线路
    circuit = QuantumCircuit(qubits)
    for i in range(len(gate_list)):
        if len(gate_list[i]) == 3: # toffoli
            circuit.ccx(gate_list[i][0],gate_list[i][1],gate_list[i][2])
        if len(gate_list[i]) == 2:  # cnot
            circuit.cx(gate_list[i][0],gate_list[i][1])
        if len(gate_list[i]) == 1:  # 单量子门 全部用x门代替
            circuit.h(gate_list[i][0])
        # 绘制线路
    # circuit.draw('mpl', scale=1, filename=img_path, initial_state=True, plot_barriers=False, justify='left', fold=-1)
    # circuit.draw('mpl',scale=1)
    print(circuit)


# 绘制DAG图
def draw_gragh(matrix, gate_list, is_global_gate_list):
    G = nx.DiGraph()
    n = len(matrix)
    point = []
    for i in range(len(gate_list)):
        gate_num = 'g' + str(i)
        point.append(gate_num)
    G.add_nodes_from(point)
    for i in range(n):
        for k in range(i + 1, n):
            if matrix[i][k] != 0:
                G.add_edge('g' + str(i), 'g' + str(k), weight=matrix[i][k])
    if Is_show_DAG == 1:
        color_map = []  # 颜色
        for j in range(len(is_global_gate_list)):
            if is_global_gate_list[j] == 0:  # 局部门
                color_map.append('green')
            if is_global_gate_list[j] == 1:  # 全局门
                color_map.append('red')
        position = nx.circular_layout(G)
        labels = nx.get_edge_attributes(G, 'weight')
        nx.draw(G, position, nodelist=point, font_color='white', node_color=color_map, with_labels=True)
        nx.draw_networkx_edge_labels(G, position, edge_labels=labels)
        plt.show()
    return G

# 输出DAG图的矩阵
def draw_matrix_by_gate_list(gate_list):
    # 量子位数
    # circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 创建对应量子位数的矩阵
    graph = Graph(len(gate_list))
    # 为矩阵添加元素
    for i in range(len(gate_list)):
        # 先找第一个量子位的下一个门
        for j in range(i + 1, len(gate_list)):
            if gate_list[i][0] in gate_list[j]:
                graph.add_edge(i, j, 'q' + str(gate_list[i][0]))
                break
        # 再找第二个量子位的下一个门
        for k in range(i + 1, len(gate_list)):
            if gate_list[i][1] in gate_list[k]:
                # 不存在其他量子位
                if graph.edge_matrix[i][k] == 0:
                    graph.add_edge(i, k, 'q' + str(gate_list[i][1]))
                    break
                # 已经存在一个依赖，再增加一个依赖
                if graph.edge_matrix[i][k] != 0:
                    graph.add_edge(i, k, graph.edge_matrix[i][k] + ',q' + str(gate_list[i][1]))
                    break
    # 矩阵输出
    # graph.print_edge_matrix()
    # 矩阵 和上面的一样
    matrix = graph.edge_matrix
    # draw_gragh(matrix)
    return matrix


'''计算每个门的标签'''
def statistics_gate_labels(gate_list,cut_list):
    statistics_gate_labels_list = [] #门标签
    partition_interval_list = [] # 分区区间 [[0,1],[2,3],[4,6]]
    c = 0
    o = 0
    while c < len(cut_list):
        # 生成区间
        partition_interval_list.append([o,o+cut_list[c]-1])
        o = o+cut_list[c]
        c += 1
    for i in range(len(gate_list)): # 每个门
        gate_labels = []
        for j in range(len(gate_list[i])): # 每个门的每个量子位
            for k in range(len(partition_interval_list)):
                if gate_list[i][j] >= partition_interval_list[k][0] and gate_list[i][j] <= partition_interval_list[k][1]:
                    gate_labels.append(k)
        statistics_gate_labels_list.append(gate_labels)
    # print(statistics_gate_labels_list)
    return statistics_gate_labels_list


'''判断gate_list中门是不是全局门'''
def is_global_gate(gate_list,cut_list):
    is_global_gate_list = []
    global_gate_num_list = []
    statistics_gate_labels_list = statistics_gate_labels(gate_list,cut_list)
    for i in range(len(statistics_gate_labels_list)):
        # 局部门
        if statistics_gate_labels_list[i][0] == statistics_gate_labels_list[i][1]:
            is_global_gate_list.append(0)
        # 全局门
        else:
            is_global_gate_list.append(1)
            global_gate_num_list.append('g'+ str(i))
    return is_global_gate_list,global_gate_num_list


# 深度优先遍历
def dfs(G, V, tags, global_gate_num_list, is_global_gate_list, statistics_gate_labels_list):
    # 记录搜索路径
    search_path = []
    # for n in G.nodes():
    #     G.nodes[n]['visited'] = False
    # print(G.nodes(data=True))
    # print('-----')
    G.nodes[V]['visited'] = True
    search_path.append(V)
    gate_label = statistics_gate_labels_list[int(V[1:])]
    travel(G, V, tags, gate_label, search_path,global_gate_num_list,is_global_gate_list,statistics_gate_labels_list)
    return search_path,G


# dfs递归实现
def travel(G, V, tags, gate_label, search_path, global_gate_num_list, is_global_gate_list, statistics_gate_labels_list):
    neighbors = nx.neighbors(G, V)
    for n in neighbors:
        # print(V + '-' + n + '-' + G.get_edge_data(V, n).get('weight'))
        # print(is_fork_in_path(G, n, search_path,is_global_gate_list))
        # 判断条件 tags in G.get_edge_data(V, n).get('weight') 边的权重相同
        # 判断条件 is_fork_in_path(G, n,search_path,is_global_gate_list) == 0 不存在依赖的全局门
        n_label = statistics_gate_labels_list[int(n[1:])]
        if not G.nodes[n]['visited'] and n in global_gate_num_list and set(n_label) == set(gate_label) and tags in G.get_edge_data(V, n).get('weight') and is_fork_in_path(G, n, search_path,is_global_gate_list) == 0:
            # 为已经访问过的节点打上标记
            G.nodes[n]['visited'] = True
            search_path.append(n)
            # 针对某个节点深入遍历搜索
            travel(G, n, tags, gate_label, search_path, global_gate_num_list, is_global_gate_list, statistics_gate_labels_list)


# 判断dfs路径中是否存在岔路
def is_fork_in_path(G, V, search_path, is_global_gate_list):
    # 两种方式实现
    fork = list(nx.ancestors(G, V))  # 不会存在V
    # T = nx.dfs_tree(G.reverse(), V)  # 存在V
    # print(list(T))
    fork_list = list(set(fork) - set(search_path))  # 差集
    #print(fork_list)
    count = []
    influ_gate = []
    for i in range(len(fork_list)):
        # print(G.nodes[fork_list[i]]['visited'])
        # 全局门 未被访问过
        if is_global_gate_list[int(fork_list[i][1:])] == 1 and not G.nodes[fork_list[i]]['visited']:
            count.append(1)
            influ_gate.append(fork_list[i])
        # 局部门
        if is_global_gate_list[int(fork_list[i][1])] == 0 :
            count.append(0)
    # print(count)
    # print(influ_gate)
    if sum(count) == 0:  # 全是局部门
        return 0
    else:
        return 1


# 生成transfer_qubit_list 初始量子位传输列表
def transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list):
    transfer_qubit_list = []
    # circuit_qubit = max([max(row) for row in gate_list]) + 1
    i = 0
    while i < len(gate_list):
        if i == len(gate_list)-1:
            transfer_qubit_list.append(gate_list[i][0])
            break
        if statistics_gate_labels_list[i][0] != statistics_gate_labels_list[i][1]:  # 是全局门
            for j in range(i + 1, len(gate_list)):
                if statistics_gate_labels_list[j][0] == statistics_gate_labels_list[j][1]:  # 下一个是局部门
                    if j == len(gate_list) - 1:
                        transfer_qubit_list.append(gate_list[i][0])
                        i = i + 1
                        break
                    else:
                        continue
                if statistics_gate_labels_list[j][0] != statistics_gate_labels_list[j][1]:  # 下一个是全局门
                    if len(list(set(gate_list[i]) & set(gate_list[j]))) == 0: # 没有相同量子位
                        if j == len(gate_list)-1:
                            transfer_qubit_list.append(gate_list[i][0])
                            i = i + 1
                            break
                        else:
                            continue
                    if len(list(set(gate_list[i]) & set(gate_list[j]))) != 0: # 有相同量子位
                        transfer_qubit_list.append(list(set(gate_list[i]) & set(gate_list[j]))[0])  # 将两个门的相同量子位加入列表
                        i = i + 1
                        break
        else: # 是局部门
            i = i + 1
    for j in range(len(transfer_qubit_list)):
        transfer_qubit_list[j] = 'q'+str(transfer_qubit_list[j])
    # print(transfer_qubit_list)
    return transfer_qubit_list

# 递归寻找合并传输队列
# G
# V传输的第一个门
# transfer_qubit_list表示所传输的量子位 = ['q0','q0','q0','q1','q0','q4','q0','q1']
# global_gate_num_list全局门序号列表 ['g0','g1','g2','g3','g4',...]
def dfs_G_to_ST(G,V,transfer_qubit_list,global_gate_num_list,ST,T_queue,is_global_gate_list,statistics_gate_labels_list):
    search_path,G = dfs(G,V,transfer_qubit_list[global_gate_num_list.index(V)],global_gate_num_list,is_global_gate_list,statistics_gate_labels_list)  # dsf(G,'g0','q0')
    # G = dfs(G,V,transfer_qubit_list[global_gate_num_list.index(V)],global_gate_num_list,is_global_gate_list,statistics_gate_labels_list)[1]  # dsf(G,'g0','q0')
    # print('search_path:'+str(search_path))
    ST = ST + search_path
    T_queue.append(search_path)
    leave_global_gate_list = list(set(global_gate_num_list)-set(ST)) # 更新全局门列表
    leave_global_gate_list = list_sort(leave_global_gate_list)
    #print(leave_global_gate_list)
    # 全局门都遍历结束
    if len(leave_global_gate_list) != 0:
        dfs_G_to_ST(G,leave_global_gate_list[0],transfer_qubit_list,global_gate_num_list,ST,T_queue,is_global_gate_list,statistics_gate_labels_list)
    return T_queue


'''list排序'''
def list_sort(str_list):
    int_list = []
    new_list = []
    for i in str_list:
        int_list.append(int(i[1:]))
    int_list.sort()
    for j in int_list:
        new_list.append('g'+str(j))
    return new_list

# 统计传输队列
def count_transfer_queue(gate_list,cut_list,transfer_qubit_list,statistics_gate_labels_list):
    # 绘制量子线路
    # new_circuit(gate_list, circuit_qubit)
    # DAG矩阵输出
    matrix = draw_matrix_by_gate_list(gate_list)
    # 判断每个门是否为全局门 全局门为 1 局部门为 0
    is_global_gate_list = is_global_gate(gate_list,cut_list)[0]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = is_global_gate(gate_list,cut_list)[1]
    # print('全局门：' + str(global_gate_num_list))
    transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)
    # 绘制DAG
    G = draw_gragh(matrix, gate_list, is_global_gate_list)
    for n in G.nodes():
        G.nodes[n]['visited'] = False
    # 变量初始化
    ST = []
    T_queue = []
    T_queue_result = dfs_G_to_ST(G, global_gate_num_list[0], transfer_qubit_list, global_gate_num_list, ST, T_queue, is_global_gate_list, statistics_gate_labels_list)
    # print(T_queue_result)
    return T_queue_result

# 考虑末尾的门后续有没有门来判断是否需要回传
def count_transfer_queue_by_consider_last(gate_list,cut_list,transfer_qubit_list,statistics_gate_labels_list):
    # 绘制量子线路
    # new_circuit(gate_list, circuit_qubit)

    # DAG矩阵输出
    matrix = draw_matrix_by_gate_list(gate_list)
    # 判断每个门是否为全局门 全局门为 1 局部门为 0
    is_global_gate_list = is_global_gate(gate_list,cut_list)[0]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = is_global_gate(gate_list,cut_list)[1]
    # print('全局门：' + str(global_gate_num_list))
    transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)
    # 绘制DAG
    G = draw_gragh(matrix, gate_list, is_global_gate_list)
    for n in G.nodes():
        G.nodes[n]['visited'] = False
    # 变量初始化
    ST = []
    T_queue = []
    T_queue_result = dfs_G_to_ST(G, global_gate_num_list[0], transfer_qubit_list, global_gate_num_list, ST, T_queue, is_global_gate_list, statistics_gate_labels_list)
    st_num = len(T_queue_result)*2
    # print(T_queue_result)
    for i in range(len(T_queue_result)):
        # 无后继节点 或者 后记节点均为执行分区的局部门
        # next_node_list = list(G.successors(T_queue_result[i][-1]))
        next_node_list = list(nx.dfs_successors(G,T_queue_result[i][-1]))
        # print(next_node_list)
        is_all_local_gate = 1
        next_node_gate_list = []
        for gate in next_node_list:
            if gate in global_gate_num_list: # 如果有一个全局门
                is_all_local_gate = 0  # 不满足均为局部门的条件
            next_node_gate_list.append(gate_list[int(gate[1:])])
        # print(next_node_gate_list)
        # 目标门的传输量子位
        # print(transfer_qubit_list)
        # print(T_queue_result[i][0])
        list_id = [element for sublist in T_queue_result for element in sublist]  # [['g1', 'g2', 'g9'], ['g10', 'g11', 'g12', 'g13', 'g14'], ['g17', 'g18']]变成一维
        target_transfer_qubit = transfer_qubit_list[list_id.index(T_queue_result[i][0])]  # 'q1'
        # print(target_transfer_qubit)
        # 均在执行分区： 目标门的传输量子位与next_node_list的量子位没有交集
        if next_node_list == [] or (is_all_local_gate == 1 and target_transfer_qubit[1:] not in list(set([element for sublist in next_node_gate_list for element in sublist]))):
            st_num -= 1/8954796859
    return len(T_queue_result)*2, st_num

# 随机扰动获取传输代价
def random_perturbations(gate_list,cut_list):
    # 量子门列表
    # gate_list = [[2, 0], [0, 2], [2, 0], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]]
    # 线路的量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 切割点 P1分区的量子位数
    # 判断每个门是否为全局门 全局门为 1 局部门为 0
    is_global_gate_list = is_global_gate(gate_list, cut_list)[0]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = is_global_gate(gate_list, cut_list)[1]
    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)

    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)

    # 统计合并传输队列
    print('######################################################################################################')
    initial_transfer_queue = count_transfer_queue(gate_list,cut_list,initial_transfer_qubit_list,statistics_gate_labels_list)
    print('初始量子位传输列表：' + str(initial_transfer_qubit_list))
    print('初始传输队列：' + str(initial_transfer_queue))
    print('初始传输代价：' + str(len(initial_transfer_queue) * 2))
    print('######################################################################################################')
    # 初始传输代价
    min_transmission_cost = len(initial_transfer_queue)*2
    # 随机扰动量子位传输列表
    for i in range(number_of_disturbance_iterations):
        transfer_qubit_list = random_perturbations_transfer_qubit_list(initial_transfer_qubit_list,gate_list,global_gate_num_list)
        print(transfer_qubit_list)
        transmission_cost = len(count_transfer_queue(gate_list,cut_list,initial_transfer_qubit_list,statistics_gate_labels_list)) * 2
        print('传输代价：' + str(transmission_cost))
        print('######################################################################################################')
        if transmission_cost < min_transmission_cost:
            min_transmission_cost = transmission_cost
    # 最低传输代价输出
    print('')
    print('***************************************************************************************************')
    print('最低传输代价：'+str(min_transmission_cost))
    print('***************************************************************************************************')
    return min_transmission_cost

# 禁忌搜索开始
# 候选者生成,生成number_of_candidate_generation个候选解,存入candidate_generation_list
def candidate_generation(gate_list, global_gate_num_list, number_of_candidate_generation, initial_transfer_qubit_list):
    candidate_generation_list = [0]*number_of_candidate_generation
    candidate_generation_list[0] = copy.deepcopy(initial_transfer_qubit_list)
    # 初始传输量子位
    for i in range(1,number_of_candidate_generation):
        transfer_qubit_list = random_perturbations_transfer_qubit_list(initial_transfer_qubit_list, gate_list, global_gate_num_list)
        candidate_generation_list[i] = copy.deepcopy(transfer_qubit_list)
    return candidate_generation_list

# 禁忌过滤 candidate_generation_list,taboo_list均为二维列表
def taboo_filter(candidate_generation_list,taboo_list):
    # print('过滤')
    # print(candidate_generation_list)
    # print(taboo_list)
    same = []
    # 在候选解列表中删除禁忌表中存在的元素
    for i in range(len(candidate_generation_list)):
        for j in range(len(taboo_list)):
            if taboo_list[j] == candidate_generation_list[i]:
                same.append(candidate_generation_list[i])
                break
    # print('same:'+str(same))
    for k in same:
        candidate_generation_list.remove(k)
    return candidate_generation_list

# 最优解选择 返回最优解对应的传输量子位列表
def optimal_solution_selection(gate_list, cut_point, candidate_generation_list, statistics_gate_labels_list):
    # 最低传输代价
    bbest_transfer_qubit_list = []
    min_transmission_cost = len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[0],statistics_gate_labels_list)) * 2
    # 最优解 最优传输量子位列表
    for i in range(len(candidate_generation_list)):
        transmission_cost = len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[i],statistics_gate_labels_list)) * 2
        if transmission_cost <= min_transmission_cost:
            # 最优解 最优传输量子位列表
            min_transmission_cost = transmission_cost
            bbest_transfer_qubit_list = copy.deepcopy(candidate_generation_list[i])
    return bbest_transfer_qubit_list


# 按传输代价对候选解排序 冒泡排序
def sort_candidate_generation(gate_list, cut_point,candidate_generation_list,statistics_gate_labels_list):
    transmission_cost_list = []
    for k in range(len(candidate_generation_list)):
        transmission_cost_list.append(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[k],statistics_gate_labels_list))*2)
    for i in range(len(candidate_generation_list)):
        for j in range(0,len(candidate_generation_list)-i-1):
            # transmission_cost_j = len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[j]))*2
            # transmission_cost_j_1 = len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[j+1])) * 2
            # if transmission_cost_j > transmission_cost_j_1:
            #     candidate_generation_list[j],candidate_generation_list[j+1] = candidate_generation_list[j+1],candidate_generation_list[j]
            if transmission_cost_list[j] > transmission_cost_list[j+1]:
                candidate_generation_list[j], candidate_generation_list[j + 1] = candidate_generation_list[j + 1], \
                                                                                 candidate_generation_list[j]
                transmission_cost_list[j],transmission_cost_list[j+1] = transmission_cost_list[j+1],transmission_cost_list[j]
    return candidate_generation_list


# 禁忌表更新 根据禁忌表大小存入最优传输代价的量子位列表，及部分传输代价较低的量子位列表
def update_taboo_list(gate_list, cut_point, best_transfer_qubit_list,taboo_list,candidate_generation_list,statistics_gate_labels_list):
    sort_candidate_generation(gate_list, cut_point, candidate_generation_list,statistics_gate_labels_list)
    # 将候选者中的最优量子位列表
    taboo_list.append(best_transfer_qubit_list)
    # 禁忌表长度超过阈值
    if len(taboo_list) > length_of_taboo_table:
        # 删除部分，留出空间给候选解
        taboo_list = sort_candidate_generation(gate_list, cut_point, taboo_list,statistics_gate_labels_list) # 排序
        # 删除 len(taboo_list)-length_of_taboo_table+1/2candidate_generation_list长度的taboo_list
        for i in range(len(taboo_list)-length_of_taboo_table+int(len(candidate_generation_list)/2)):
            taboo_list.pop()
        # taboo_list补全到length_of_taboo_table个元素
        ik = 1
        while len(taboo_list) < length_of_taboo_table:
            taboo_list.append(candidate_generation_list[ik])
            ik += 1
    # 禁忌表长度小于阈值
    if len(taboo_list) < length_of_taboo_table:
        j = 1
        while len(taboo_list) < length_of_taboo_table:
            taboo_list.append(candidate_generation_list[j])
            j += 1
            if j >= len(candidate_generation_list):
                break
    # 禁忌表去重
    new_taboo_list = []
    for k in taboo_list:
        if k not in new_taboo_list:
            new_taboo_list.append(k)
    return new_taboo_list


# 禁忌搜索获取传输代价
def taboo_search(gate_list,cut_list):
    # 线路的量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 判断每个门是否为全局门 全局门为 1 局部门为 0
    is_global_gate_list = is_global_gate(gate_list, cut_list)[0]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = is_global_gate(gate_list, cut_list)[1]

    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)
    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)

    # 统计合并传输队列
    print('######################################################################################################')
    initial_transfer_queue = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list,
                                                  statistics_gate_labels_list)
    print('初始量子位传输列表：' + str(initial_transfer_qubit_list))
    print('全局门：'+str(global_gate_num_list))
    print('初始传输队列：' + str(initial_transfer_queue))
    print('初始传输代价：' + str(len(initial_transfer_queue) * 2))
    print('######################################################################################################')

    taboo_list = []  # 禁忌表
    # 开始迭代
    for i in range(number_of_iterations):
        print('迭代第' + str(i) + '次,最优传输代价：',end=' ')
        # 生成候选解
        candidate_generation_list = candidate_generation(gate_list, global_gate_num_list, number_of_candidate_generation, initial_transfer_qubit_list)
        # 禁忌过滤
        candidate_generation_list = taboo_filter(candidate_generation_list, taboo_list)
        # candidate_generation_list = sort_candidate_generation(candidate_generation_list)
        # print(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[0]))*2)
        # print('禁忌过滤后的候选解' + str(candidate_generation_list))
        transmission_cost_list = []
        for j in range(len(candidate_generation_list)):
            transmission_cost_j = len(count_transfer_queue(gate_list, cut_list, candidate_generation_list[j],statistics_gate_labels_list)) * 2
            print(transmission_cost_j, end=' ')
            transmission_cost_list.append(transmission_cost_j)
        # 最优解生成
        best_transfer_qubit_list = optimal_solution_selection(gate_list, cut_list,candidate_generation_list,statistics_gate_labels_list)
        # best_transfer_qubit_list = candidate_generation_list[transmission_cost_list.index(min(transmission_cost_list))]
        print(' ', end=' ')
        print(len(count_transfer_queue(gate_list, cut_list, best_transfer_qubit_list, statistics_gate_labels_list)) * 2)
        # 禁忌表更新
        # print('候选解中最优解' + str(best_transfer_qubit_list) + ' 传输代价：' + str(len(count_transfer_queue(gate_list, cut_point, best_transfer_qubit_list)) * 2))
        taboo_list = update_taboo_list(gate_list, cut_list,best_transfer_qubit_list, taboo_list, candidate_generation_list,statistics_gate_labels_list)

        # print('禁忌表：' + str(taboo_list))
        # print('#############################################################')
    # 禁忌表最优解输出
    best_transfer_qubit_list_taboo = optimal_solution_selection(gate_list, cut_list, taboo_list,statistics_gate_labels_list)
    min_T_queue_list = count_transfer_queue(gate_list, cut_list, best_transfer_qubit_list_taboo,statistics_gate_labels_list)
    min_transmission_cost = len(min_T_queue_list) * 2
    print('######################################################################################################')
    print('最低传输代价：' + str(min_transmission_cost))
    print('最优量子位传输列表：' + str(best_transfer_qubit_list_taboo))
    print('最优传输队列：' + str(min_T_queue_list))
    print('######################################################################################################')
    return min_transmission_cost


# 禁忌搜索获取传输代价 无输出版
def taboo_search2(gate_list,cut_list):
    # 线路的量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 判断每个门是否为全局门 全局门为 1 局部门为 0
    is_global_gate_list = is_global_gate(gate_list, cut_list)[0]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = is_global_gate(gate_list, cut_list)[1]
    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)
    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)

    taboo_list = []  # 禁忌表
    # 开始迭代
    for i in range(number_of_iterations):
        # 生成候选解
        candidate_generation_list = candidate_generation(gate_list, global_gate_num_list, number_of_candidate_generation, initial_transfer_qubit_list)
        # 禁忌过滤
        candidate_generation_list = taboo_filter(candidate_generation_list, taboo_list)
        # candidate_generation_list = sort_candidate_generation(candidate_generation_list)
        # print(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[0]))*2)
        # print('禁忌过滤后的候选解' + str(candidate_generation_list))
        # for j in range(len(candidate_generation_list)):
        #     print(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[j])) * 2, end=' ')
        # 最优解生成
        best_transfer_qubit_list = optimal_solution_selection(gate_list, cut_list,candidate_generation_list,statistics_gate_labels_list)
        # 禁忌表更新
        # print('候选解中最优解' + str(best_transfer_qubit_list) + ' 传输代价：' + str(len(count_transfer_queue(gate_list, cut_point, best_transfer_qubit_list)) * 2))
        taboo_list = update_taboo_list(gate_list, cut_list,best_transfer_qubit_list, taboo_list, candidate_generation_list,statistics_gate_labels_list)
        # print('禁忌表：' + str(taboo_list))
    # 禁忌表最优解输出
    best_transfer_qubit_list_taboo = optimal_solution_selection(gate_list, cut_list, taboo_list,statistics_gate_labels_list)
    min_T_queue_list = count_transfer_queue(gate_list, cut_list, best_transfer_qubit_list_taboo,statistics_gate_labels_list)
    min_transmission_cost = len(min_T_queue_list) * 2
    print(str(gate_list)+' '+str(min_T_queue_list)+' 传输代价：'+str(min_transmission_cost))
    return min_transmission_cost


'''初始解直接搜索'''
def direct_calculation_of_tc(gate_list,cut_list):
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
#    global_gate_num_list = is_global_gate(gate_list, cut_list)[1]

    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)
    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)

    # 统计合并传输队列
    #print('######################################################################################################')
    initial_transfer_queue = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list,
                                                  statistics_gate_labels_list)
    # print('初始量子位传输列表：' + str(initial_transfer_qubit_list))
    # print('全局门：' + str(global_gate_num_list))
    # print('初始传输队列：' + str(initial_transfer_queue))
    print('初始传输代价：' + str(len(initial_transfer_queue) * 2),end=' ')
    # print('######################################################################################################')

    return len(initial_transfer_queue) * 2



'''初始解直接搜索 通过队列'''
def direct_calculation_of_tc_queue(gate_list,cut_list):
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
#    global_gate_num_list = is_global_gate(gate_list, cut_list)[1]
#     print(gate_list)
    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)
    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, statistics_gate_labels_list)
    # 全局门列表['g2', 'g6', 'g9']
    global_gate_list = new_global_gate_list(statistics_gate_labels_list)
    # print(global_gate_list)
    # 统计合并传输队列
    #print('######################################################################################################')
    initial_transfer_queue = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list,
                                                  statistics_gate_labels_list)
    # print('初始量子位传输列表：' + str(initial_transfer_qubit_list))
    # print('全局门：' + str(global_gate_num_list))
    # print('初始传输队列：' + str(initial_transfer_queue))
    st = len(initial_transfer_queue) * 2
    print('初始传输代价：' + str(st),end=' ')
    # print('######################################################################################################')

    # 如果两个队列间隔的一些量子门与传输量子位无相同量子位，且量子位传输列表位于同一个分区，则回传和正传可省略
    for i in range(len(initial_transfer_queue)-1):

        # 量子位传输列表位于同一个分区
        gate_i = initial_transfer_queue[i][0]  # g2
        gate_num_i = global_gate_list.index(gate_i) # ['g2', 'g6', 'g9'] -> 0
        transfer_qubit_i = int(get_data(initial_transfer_qubit_list[gate_num_i])[0]) #
        gate_j = initial_transfer_queue[i+1][0]  # g6
        gate_num_j = global_gate_list.index(gate_j)
        transfer_qubit_j = int(get_data(initial_transfer_qubit_list[gate_num_j])[0])

        gate_i_next_name = initial_transfer_queue[i][-1]  # ['g2', 'g6', 'g9'] -> g9
        gate_i_next = int(gate_i_next_name[1:])+1  # [[2, 0], [0, 2], [2, 5], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]] -> 3
        gate_j_front = int(gate_j[1:])  # -> 6
        middle_gate_list = gate_list[gate_i_next:gate_j_front]  # 左闭右开，所以 gate_j_front 不用 +1

        # 两个队列间隔的一些量子门与传输量子位无相同量子位   起始点：第一个队列中第一个门 结束点：第二个队列中第一个门
        if judge_is_same_partation(transfer_qubit_i,transfer_qubit_j,cut_list) == 1 and judge_middle_gates_is_affect_transfer_qubit(middle_gate_list,transfer_qubit_i)==0:
            # print(gate_i)
            # print(gate_j)
            st -= 1
    print('优化后传输代价：' + str(st),end=' ')
    return st



# 全局门列表
def new_global_gate_list(statistics_gate_labels_list):
    global_gate_list = []
    for i in range(len(statistics_gate_labels_list)):
        if statistics_gate_labels_list[i][0] != statistics_gate_labels_list[i][1]:
            global_gate_list.append('g'+str(i))
    return global_gate_list


# 判断两个量子位是否位于同分区
def judge_is_same_partation(transfer_qubit_i,transfer_qubit_j,cut_list):
    partition_interval_list = []  # 分区区间 [[0,1],[2,3],[4,6]]
    c = 0
    o = 0
    while c < len(cut_list):
        # 生成区间
        partition_interval_list.append([o, o + cut_list[c] - 1])
        o = o + cut_list[c]
        c += 1
    for interval in partition_interval_list:
        if interval[0] <= transfer_qubit_i <= interval[1] and interval[0] <= transfer_qubit_j <= interval[1]:
            return 1
    return 0

# 判断一些量子门是否会影响到传输量子位
def judge_middle_gates_is_affect_transfer_qubit(middle_gate_list,transfer_qubit):
    for gate in middle_gate_list:
        if transfer_qubit in gate:
            return 1
    return 0



'''初始解直接搜索-分布式方法 N：线路分分成几份'''
def direct_calculation_of_tc_distributed(gate_list,cut_list,N):
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    gate_list_list = split_list(gate_list, N)

    # 统计合并传输队列
    #print('######################################################################################################')
    # initial_transfer_queue = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list,statistics_gate_labels_list)
    initial_transfer_queue_list = []
    for i in range(N):
        # 全局门标签
        statistics_gate_labels_list = statistics_gate_labels(gate_list_list[i], cut_list)
        # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
        initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list_list[i], statistics_gate_labels_list)
        try:
            initial_transfer_queue = count_transfer_queue(gate_list_list[i], cut_list, initial_transfer_qubit_list,
                                                          statistics_gate_labels_list)
            initial_transfer_queue_list.append(initial_transfer_queue)
            # print(len(initial_transfer_queue)*2)
        except Exception as e:
            # print(gate_list_list[i])
            pass

    st = 0
    for j in range(len(initial_transfer_queue_list)):
        st += len(initial_transfer_queue_list[j])
    print('初始传输代价：' + str(st * 2))

    return st * 2


# 将list分割成N份
def split_list(lst, n):
    # 计算每份的大小
    avg = len(lst) // n
    remainder = len(lst) % n

    # 初始化结果列表
    result = []

    # 分割列表
    start = 0
    for i in range(n):
        end = start + avg + (1 if i < remainder else 0)
        result.append(lst[start:end])
        start = end

    return result


if __name__ == '__main__':
    # 门列表
    # gate_list = [[4, 3], [7, 4], [3, 7], [3, 4], [7, 4], [3, 7], [4, 3], [7, 6], [1, 7], [6, 1], [6, 7], [1, 7], [6, 1], [7, 6], [4, 3], [2, 4], [3, 2], [3, 4], [2, 4], [3, 2], [4, 3], [4, 2], [2, 0], [1, 2], [0, 1], [0, 2], [1, 2], [0, 1], [2, 0], [1, 0], [2, 1], [0, 2], [0, 1], [2, 1], [0, 2], [1, 0], [3, 5], [4, 5], [6, 8], [7, 8], [5, 6], [8, 5], [6, 8], [6, 5], [8, 5], [6, 8], [5, 6], [7, 6], [8, 7], [6, 8], [6, 7], [8, 7], [6, 8], [7, 6], [6, 9], [7, 6], [9, 7], [6, 9], [6, 7], [9, 7], [6, 9], [7, 6], [7, 9], [9, 0], [8, 9], [0, 8], [0, 9], [8, 9], [0, 8], [9, 0], [8, 0], [9, 8], [0, 9], [0, 8], [9, 8], [0, 9], [8, 0]]
    # gate_list = [[2, 0], [0, 2], [2, 5], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]]
    gate_list = [[0,3],[5,3],[3,5],[0,3],[0,5],[3,5],[5,3],[1,5],[4,5],[5,4],[1,5],[1,4],[5,4],[4,5],[1,3],[2,3],[3,2],[1,3],[1,2],[3,2],[2,3]]
    draw_new_circuit(gate_list)
    # random_perturbations(gate_list, cut_list)

    # taboo_search(gate_list, cut_list)
    cut_list = [2,2,2]
    direct_calculation_of_tc_queue(gate_list,cut_list)


