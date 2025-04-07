import datetime

from Utils.division_by_irregularity import direct_calculation_of_tc_look_ahead
from Utils.generate_partitions import generate_line
from Utils.min_global_gate_num import k_change_gate_list_by_line_sequence, group_nodes_by_partition, \
    Initial_line_sequence, random_line_by_qubit
from Utils.read_qasm import list_str_to_int, converter_circ_from_qasm
from Utils.transmission_cost_calculation_more import is_global_gate, direct_calculation_of_tc

import sys

from partition.matrix_to_pujulei import gate_list_to_graph, spectral_min_cut

sys.setrecursionlimit(50000)  # 设置更高的递归深度限制



# 根据划分方法输出划分方案的labels
def count_st_by_partition(partition_method,k,gate_list):
    if partition_method == 'spectral_clustering':
        graph = gate_list_to_graph(gate_list)
        partition_labels = spectral_min_cut(graph,k)
        st = k_count_st_num_ahead(partition_labels, gate_list)

    return st



# 按特定划分算法下的传输代价
# partition_labels: [0, 2, 2, 1, 1, 0, 0, 1, 1, 1] 每个节点所属的分区编号
def k_count_st_num_ahead(partition_labels, gate_list):

    node_labels = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z','0','1','2','3','4','5','6','7','8','9']

    # cut_list: [3, 5, 2] 每个分区大小 ; combined_line_sequence: "AFGDEHIJBC" 线序
    cut_list, combined_line_sequence = group_nodes_by_partition(partition_labels, node_labels[:len(partition_labels)])

    initial_line_sequence = Initial_line_sequence(gate_list)
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list, combined_line_sequence)
    st = direct_calculation_of_tc_look_ahead(new_gate_list, cut_list)

    return st


# 随机划分
def k_count_min_st_num_ahead_by_random_partition(gate_list, cut_list, random_num):

    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    line_sequence_list = random_line_by_qubit(random_num,circuit_qubit)
    lenl = len(line_sequence_list)
    print(lenl)
    mini = len(gate_list) * 2
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i], end=' ')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence_list[i])
        st = direct_calculation_of_tc_look_ahead(new_gate_list, cut_list)
        ST.append(st)
        print('当前最小：', mini)
        if i % 100 == 0: # 100次统计一次
            mini = min(ST)

    print('最小：', min(ST))
    # min_st_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list,
    #                                                        line_sequence_list[ST.index(min(ST))])
    # 绘制电路
    # draw_color_circuit(min_st_gate_list)
    return ST


# 会遍历所有划分情况
def k_count_min_st_num_ahead_by_iteration(gate_list,cut_list):
    line_sequence_list = generate_line(cut_list)
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = [] # 传输代价
    lenl = len(line_sequence_list)

    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end=' ')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])

        global_gate_num_list = is_global_gate(new_gate_list, cut_list)[1]
        print('全局门：', end='')
        print(len(global_gate_num_list), end=' ')

        print('当前线序下前瞻的', end='')

        st0 = direct_calculation_of_tc(new_gate_list, cut_list) # 无前瞻
        st_ahead = direct_calculation_of_tc_look_ahead(new_gate_list, cut_list)  # 前瞻
        print()
        ST.append(st_ahead)

    print('最小', min(ST))

    return ST



if __name__ == '__main__':

    start_time = datetime.datetime.now()

    input_filename = r'./qasm/qasm_czl/mini_alu_305.qasm'
    # gate_list = real_to_cnot_list(input_filename)
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    print(gate_list)

    cut_list = [5,5]
    line_sequence_list = generate_line(cut_list)
    print(line_sequence_list)

    circuit_qubit = max([max(row) for row in gate_list]) + 1

    # st = k_count_st_num_ahead(partition_labels, gate_list)

    st = k_count_min_st_num_ahead_by_iteration(gate_list, cut_list)

    st = k_count_min_st_num_ahead_by_random_partition(gate_list, cut_list, 100)

    # print(st)

    end_time = datetime.datetime.now()
    print(end_time-start_time)

