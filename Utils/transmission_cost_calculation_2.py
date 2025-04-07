'''
功能：将线路转换为依赖图的矩阵形式，并绘制依赖图。根据依赖图计算分布式线路的传输代价。
注意：此功能只适用于两个分区的分布式线路的传输代价计算。
开始时间：2023.2.7
作者：陈新宇
版本: 1.0
完成时间：2023.2.22
'''

import copy
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from qiskit import QuantumRegister, ClassicalRegister, AncillaRegister, QuantumCircuit
import random

# 随机扰动迭代次数
number_of_disturbance_iterations = 100

# 禁忌搜索参数设置
number_of_iterations = 10  # 迭代次数
number_of_times_not_updated = 50  # 最优解不更新次数 达到后直接返回禁忌表中最优解
length_of_taboo_table = 10  # 禁忌表长度
number_of_candidate_generation = 10  # 候选解数量

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


# 绘制DAG图
def draw_gragh(matrix, gate_list, is_global_gate_list,is_show_graph):
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
    if is_show_graph == 1:
        plt.show()
    return G


# 输出DAG图的矩阵
def draw_matrix_by_gate_list(gate_list):
    # 量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
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

# 绘制线路
def new_circuit(gate_list,circuit_qubit):
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

'''判断一个门是否为全局门'''
def judge_is_global_gate(gate,cut_point,circuit_qubit):
    cut_list_one = []
    cut_list_two = []
    # 第一个分区的所有量子位集合
    for i in range(cut_point):
        cut_list_one.append(i)
    # 第二个分区的所有量子位集合
    for i in range(cut_point, circuit_qubit):
        cut_list_two.append(i)
    res1 = list(set(gate) & set(cut_list_one))
    res2 = list(set(gate) & set(cut_list_two))
    if (len(res1) == len(set(gate)) and (len(res2) == 0)) or (
            len(res2) == len(set(gate)) and (len(res1) == 0)):
        is_global_gate = 0
    else:
        is_global_gate = 1
    return is_global_gate


'''2分区中计算全局门个数 （门列表，切线位置，量子位数）'''
def count_global_gate_num_2(gate_list, cut_point):
    global_gate_list = []
    local_gate_list = []
    count_num = 0
    cut_list_one = []
    cut_list_two = []
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    global_gate_num_list = []  # 全局门序号列表 ['g0','g3','g4',...]
    # 第一个分区的所有量子位集合
    for i in range(cut_point):
        cut_list_one.append(i)
    # print(cut_list_one)
    # 第二个分区的所有量子位集合
    for i in range(cut_point, circuit_qubit):
        cut_list_two.append(i)
    # print(cut_list_two)
    is_global_gate_list = [0] * len(gate_list)
    # 门列表与分区的交集为门列表本身，并且与另一个分区没有交集，是局部门，否则为全局门
    for i in range(len(gate_list)):
        res1 = list(set(gate_list[i]) & set(cut_list_one))
        res2 = list(set(gate_list[i]) & set(cut_list_two))
        # 局部门
        if (len(res1) == len(set(gate_list[i])) and (len(res2) == 0)) or (
                len(res2) == len(set(gate_list[i])) and (len(res1) == 0)):
            local_gate_list.append(gate_list[i])
        # 全局门
        else:
            global_gate_list.append(gate_list[i])
            global_gate_num_list.append('g'+str(i))
            is_global_gate_list[i] = 1
            count_num = count_num + 1
    # print(global_gate_list)
    # print(local_gate_list)
    return count_num, global_gate_list, global_gate_num_list, is_global_gate_list


# 深度优先遍历
def dfs(G, V, tags, global_gate_num_list,is_global_gate_list):
    # 记录搜索路径
    search_path = []
    # for n in G.nodes():
    #     G.nodes[n]['visited'] = False
    # print(G.nodes(data=True))
    # print('-----')
    G.nodes[V]['visited'] = True
    search_path.append(V)
    travel(G, V, tags, search_path,global_gate_num_list,is_global_gate_list)
    return search_path,G


# dfs递归实现
def travel(G, V, tags, search_path,global_gate_num_list,is_global_gate_list):
    neighbors = nx.neighbors(G, V)
    for n in neighbors:
        # print(V + '-' + n + '-' + G.get_edge_data(V, n).get('weight'))
        # print(is_fork_in_path(G, n, search_path,is_global_gate_list))
        # 判断条件 tags in G.get_edge_data(V, n).get('weight') 边的权重相同
        # 判断条件 is_fork_in_path(G, n,search_path,is_global_gate_list) == 0 不存在依赖的全局门
        if not G.nodes[n]['visited'] and n in global_gate_num_list and tags in G.get_edge_data(V, n).get('weight') \
                and is_fork_in_path(G, n, search_path,is_global_gate_list) == 0:
            # 为已经访问过的节点打上标记
            G.nodes[n]['visited'] = True
            search_path.append(n)
            # 针对某个节点深入遍历搜索
            travel(G, n, tags, search_path, global_gate_num_list,is_global_gate_list)

# 递归寻找合并传输队列
# G,
# V传输的第一个门,
# transfer_qubit_list表示所传输的量子位 = ['q0','q0','q0','q1','q0','q4','q0','q1']
# global_gate_num_list全局门序号列表 ['g0','g1','g2','g3','g4',...]
def dfs_G_to_ST(G,V,transfer_qubit_list,global_gate_num_list,ST,T_queue,is_global_gate_list):
    search_path = dfs(G,V,transfer_qubit_list[global_gate_num_list.index(V)],global_gate_num_list,is_global_gate_list)[0]  # dsf(G,'g0','q0')
    G = dfs(G,V,transfer_qubit_list[global_gate_num_list.index(V)],global_gate_num_list,is_global_gate_list)[1]  # dsf(G,'g0','q0')
    # print('search_path:'+str(search_path))
    ST = ST + search_path
    T_queue.append(search_path)
    leave_global_gate_list = list(set(global_gate_num_list)-set(ST)) # 更新全局门列表
    leave_global_gate_list.sort() # 排序
    # print(leave_global_gate_list)
    # 全局门都遍历结束
    if len(leave_global_gate_list) != 0:
        dfs_G_to_ST(G,leave_global_gate_list[0],transfer_qubit_list,global_gate_num_list,ST,T_queue,is_global_gate_list)
    return T_queue


# 根据门获取transfer_qubit_list中传输量子位 gate='g3'
def transfer_qubit_by_gate(transfer_qubit_list,gate,global_gate_num_list):
    gate_num = gate[1:] #提取数字 'g3'->3
    gate_location_num = global_gate_num_list.index(gate) #索引位置
    return transfer_qubit_list[gate_location_num]


# 判断dfs路径中是否存在岔路
def is_fork_in_path(G, V, search_path, is_global_gate_list):
    # 两种方式实现
    fork = list(nx.ancestors(G, V))  # 不会存在V
    # T = nx.dfs_tree(G.reverse(), V)  # 存在V
    # print(list(T))
    fork_list = list(set(fork) - set(search_path))  # 差集
    # print(fork_list)
    count = []
    for i in range(len(fork_list)):
        # print(G.nodes[fork_list[i]]['visited'])
        # 全局门
        if is_global_gate_list[int(fork_list[i][1])] == 1 and not G.nodes[fork_list[i]]['visited']:
            count.append(1)
        # 局部门
        if is_global_gate_list[int(fork_list[i][1])] == 0 :
            count.append(0)
    #print(sum(count))
    if sum(count) == 0:  # 全是局部门
        return 0
    else:
        return 1


# 生成transfer_qubit_list
def transfer_qubit_list_by_gate_list(gate_list, cut_point):
    transfer_qubit_list = []
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    i = 0
    while i < len(gate_list):
        if i == len(gate_list)-1:
            transfer_qubit_list.append(gate_list[i][0])
            break
        if judge_is_global_gate(gate_list[i], cut_point,circuit_qubit) == 1:  # 是全局门
            for j in range(i + 1, len(gate_list)):
                if judge_is_global_gate(gate_list[j], cut_point,circuit_qubit) == 1:  # 是全局门
                    if len(list(set(gate_list[i]) & set(gate_list[j]))) != 0: # 有相同量子位
                        transfer_qubit_list.append(list(set(gate_list[i]) & set(gate_list[j]))[0])  # 将两个门的相同量子位加入列表
                        i = i + 1
                        break
                    else: # 没有相同量子位
                        transfer_qubit_list.append(gate_list[i][0])  # 将两个门的相同量子位加入列表
                        i = i + 1
                        break
                if judge_is_global_gate(gate_list[j], cut_point, circuit_qubit) == 0:  # 是局部门
                    transfer_qubit_list.append(gate_list[i][0]) # 将两个门的相同量子位加入列表
                    i = i + 1
                    break
        else: # 是局部门
            i = i + 1
    for j in range(len(transfer_qubit_list)):
        transfer_qubit_list[j] = 'q'+str(transfer_qubit_list[j])
    # print(transfer_qubit_list)
    return transfer_qubit_list

#随机扰动 transfer_qubit_list
def random_perturbations_transfer_qubit_list(transfer_qubit_list,gate_list,global_gate_num_list):
    for i in range(len(global_gate_num_list)):
        random_num = random.randint(0, 1)
        transfer_qubit_list[i] = 'q' + str(gate_list[int(global_gate_num_list[i][1:])][random_num])
    return transfer_qubit_list

# 统计传输队列
def count_transfer_queue(gate_list,cut_point,transfer_qubit_list):
    # 绘制量子线路
    # new_circuit(gate_list, circuit_qubit)
    # DAG矩阵输出
    matrix = draw_matrix_by_gate_list(gate_list)
    # 判断每个门是否为全局门 全局门为1 局部门为0
    is_global_gate_list = count_global_gate_num_2(gate_list, cut_point)[-1]
    global_gate_num_list = count_global_gate_num_2(gate_list, cut_point)[2]
    # print('全局门：' + str(global_gate_num_list))
    transfer_qubit_list_by_gate_list(gate_list, cut_point)
    # 绘制DAG
    G = draw_gragh(matrix, gate_list, is_global_gate_list,0)
    for n in G.nodes():
        G.nodes[n]['visited'] = False
    # 变量初始化
    ST = []
    T_queue = []
    T_queue_result = dfs_G_to_ST(G, global_gate_num_list[0], transfer_qubit_list, global_gate_num_list, ST, T_queue, is_global_gate_list)
    # print(T_queue_result)
    return T_queue_result


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
def optimal_solution_selection(gate_list, cut_point,candidate_generation_list):
    # 最低传输代价
    bbest_transfer_qubit_list = []
    min_transmission_cost = len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[0])) * 2
    # 最优解 最优传输量子位列表
    for i in range(len(candidate_generation_list)):
        transmission_cost = len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[i])) * 2
        if transmission_cost <= min_transmission_cost:
            # 最优解 最优传输量子位列表
            min_transmission_cost = transmission_cost
            bbest_transfer_qubit_list = copy.deepcopy(candidate_generation_list[i])
    return bbest_transfer_qubit_list

# 按传输代价对候选解排序 冒泡排序
def sort_candidate_generation(gate_list, cut_point,candidate_generation_list):
    transmission_cost_list = []
    for k in range(len(candidate_generation_list)):
        transmission_cost_list.append(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[k]))*2)
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
def update_taboo_list(gate_list, cut_point, best_transfer_qubit_list,taboo_list,candidate_generation_list):
    sort_candidate_generation(gate_list, cut_point, candidate_generation_list)
    # 将候选者中的最优量子位列表
    taboo_list.append(best_transfer_qubit_list)
    # 禁忌表长度超过阈值
    if len(taboo_list) > length_of_taboo_table:
        # 删除部分，留出空间给候选解
        taboo_list = sort_candidate_generation(gate_list, cut_point, taboo_list) # 排序
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

'''更新量子位传输列表'''
def update_transfer_qubit_list(transfer_qubit_list,transfer_queue):
    count = 0
    for i in range(len(transfer_queue)):
        count += len(transfer_queue[i])
        if len(transfer_queue[i]) == 1:
            continue
        if len(transfer_queue[i]) > 1:
            for j in range(count-len(transfer_queue[i])+1,count):
                transfer_qubit_list[j] = copy.deepcopy(transfer_qubit_list[count-len(transfer_queue[i])])
    return transfer_qubit_list


# 随机扰动获取传输代价
def random_perturbations(gate_list,cut_point):
    # 量子门列表
    # gate_list = [[2, 0], [0, 2], [2, 0], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]]
    # 线路的量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 切割点 P1分区的量子位数
    # cut_point = 2
    # 判断每个门是否为全局门 全局门为 1 局部门为 0
    is_global_gate_list = count_global_gate_num_2(gate_list, cut_point)[-1]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = count_global_gate_num_2(gate_list, cut_point)[2]
    # 量子位传输列表 ['q0','q0','q2','q0','q4','q1','q0','q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, cut_point)

    # 统计合并传输队列
    print('######################################################################################################')
    initial_transfer_queue = count_transfer_queue(gate_list, cut_point, initial_transfer_qubit_list)
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
        transmission_cost = len(count_transfer_queue(gate_list, cut_point, transfer_qubit_list))*2
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


# 禁忌搜索获取传输代价
def taboo_search(gate_list,cut_point):
    # 线路的量子位数
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    # 判断每个门是否为全局门 全局门为1 局部门为0
    is_global_gate_list = count_global_gate_num_2(gate_list, cut_point)[-1]
    # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
    global_gate_num_list = count_global_gate_num_2(gate_list, cut_point)[2]
    # 量子位传输列表 ['q0','q0','q2','q0','q4','q1','q0','q4']
    initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, cut_point)

    # 统计合并传输队列
    print('######################################################################################################')
    initial_transfer_queue = count_transfer_queue(gate_list, cut_point, initial_transfer_qubit_list)
    print('初始量子位传输列表：' + str(initial_transfer_qubit_list))
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
        # for j in range(len(candidate_generation_list)):
        #     print(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[j])) * 2, end=' ')
        # 最优解生成
        best_transfer_qubit_list = optimal_solution_selection(gate_list, cut_point,candidate_generation_list)
        # 禁忌表更新
        # print('候选解中最优解' + str(best_transfer_qubit_list) + ' 传输代价：' + str(len(count_transfer_queue(gate_list, cut_point, best_transfer_qubit_list)) * 2))
        taboo_list = update_taboo_list(gate_list, cut_point,best_transfer_qubit_list, taboo_list, candidate_generation_list)
        print(len(count_transfer_queue(gate_list, cut_point, optimal_solution_selection(gate_list, cut_point, taboo_list))) * 2)
        # print('禁忌表：' + str(taboo_list))
        # print('#############################################################')
    # 禁忌表最优解输出
    best_transfer_qubit_list_taboo = optimal_solution_selection(gate_list, cut_point, taboo_list)
    min_T_queue_list = count_transfer_queue(gate_list, cut_point, best_transfer_qubit_list_taboo)
    min_transmission_cost = len(min_T_queue_list) * 2
    print('######################################################################################################')
    print('最低传输代价：' + str(min_transmission_cost))
    print('最优量子位传输列表：' + str(update_transfer_qubit_list(best_transfer_qubit_list_taboo,min_T_queue_list)))
    print('最优传输队列：' + str(min_T_queue_list))
    print('######################################################################################################')
    return min_transmission_cost


if __name__ == '__main__':
    gate_list = [[2, 0], [0, 2], [2, 0], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]]
    cut_point = 2
    random_perturbations(gate_list, cut_point)

# # 禁忌搜索整个流程
# if __name__ == '__main__':
#     # gate_list = [[0, 3], [2, 3], [3, 2], [0, 3], [0, 2], [3, 2], [2, 3], [1, 2], [4, 2], [2, 4], [1, 2], [1, 4], [2, 4],[4, 2], [1, 3], [5, 3], [3, 5], [1, 3], [1, 5], [3, 5], [5, 3]]
#     # 量子门列表
#     gate_list = [[2, 0], [0, 2], [2, 0], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]]
#     # 线路的量子位数
#     circuit_qubit = max([max(row) for row in gate_list]) + 1
#     # 切割点 P1分区的量子位数
#     cut_point = 2
#     # 判断每个门是否为全局门 全局门为1 局部门为0
#     is_global_gate_list = count_global_gate_num_2(gate_list, cut_point)[-1]
#     # 全局门列表 ['g0', 'g1', 'g2', 'g4', 'g5', 'g6', 'g8', 'g9']
#     global_gate_num_list = count_global_gate_num_2(gate_list, cut_point)[2]
#     # 量子位传输列表 ['q0','q0','q2','q0','q4','q1','q0','q4']
#     initial_transfer_qubit_list = transfer_qubit_list_by_gate_list(gate_list, cut_point)
#
#     # 统计合并传输队列
#     print('######################################################################################################')
#     initial_transfer_queue = count_transfer_queue(gate_list, cut_point, initial_transfer_qubit_list)
#     print('初始量子位传输列表：' + str(initial_transfer_qubit_list))
#     print('初始传输队列：'+ str(initial_transfer_queue))
#     print('初始传输代价：'+ str(len(initial_transfer_queue)*2))
#     print('######################################################################################################')
#     # # 初始传输代价
#     # min_transmission_cost = len(initial_transfer_queue)*2
#     # # 随机扰动量子位传输列表
#     # for i in range(number_of_iterations):
#     #     transfer_qubit_list = random_perturbations_transfer_qubit_list(transfer_qubit_list,gate_list,global_gate_num_list)
#     #     transmission_cost = len(count_transfer_queue(gate_list, cut_point, transfer_qubit_list))*2
#     #     print('传输代价：' + str(transmission_cost))
#     #     print('######################################################################################################')
#     #     if transmission_cost < min_transmission_cost:
#     #         min_transmission_cost = transmission_cost
#     # # 最低传输代价输出
#     # print('')
#     # print('***************************************************************************************************')
#     # print('最低传输代价：'+str(min_transmission_cost))
#     # print('***************************************************************************************************')
#
#     # 禁忌搜索参数设置
#     number_of_iterations = 100  # 迭代次数
#     number_of_times_not_updated = 50  # 最优解不更新次数 达到后直接返回禁忌表中最优解
#     length_of_taboo_table = 10  # 禁忌表长度
#     number_of_candidate_generation = 10 # 候选解数量
#     taboo_list = [] # 禁忌表
#     # 开始迭代
#     for i in range(number_of_iterations):
#         print('迭代：第'+str(i)+'次')
#         # 生成候选解
#         candidate_generation_list = candidate_generation(number_of_candidate_generation, initial_transfer_qubit_list)
#         # 禁忌过滤
#         candidate_generation_list = taboo_filter(candidate_generation_list, taboo_list)
#         # candidate_generation_list = sort_candidate_generation(candidate_generation_list)
#         # print(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[0]))*2)
#         print('禁忌过滤后的候选解' + str(candidate_generation_list))
#         for j in range(len(candidate_generation_list)):
#             print(len(count_transfer_queue(gate_list, cut_point, candidate_generation_list[j]))*2,end=' ')
#
#         # 最优解生成
#         best_transfer_qubit_list = optimal_solution_selection(candidate_generation_list)
#         # 禁忌表更新
#         print('候选解中最优解'+str(best_transfer_qubit_list)+' 传输代价：'+str(len(count_transfer_queue(gate_list, cut_point, best_transfer_qubit_list))*2))
#         taboo_list = update_taboo_list(best_transfer_qubit_list, taboo_list, candidate_generation_list)
#         print('禁忌表：'+str(taboo_list))
#         print('#############################################################')
#     # 禁忌表最优解输出
#     best_transfer_qubit_list_taboo = optimal_solution_selection(taboo_list)
#     print('最低传输代价'+ str(len(count_transfer_queue(gate_list, cut_point, best_transfer_qubit_list_taboo))*2))








