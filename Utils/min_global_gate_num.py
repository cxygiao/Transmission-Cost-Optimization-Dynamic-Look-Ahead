'''
最小化全局门数量
陈新宇
2023.7.9
'''
import datetime
import re
import time
from multiprocessing.pool import Pool
import itertools
import random
import Utils.transmission_cost_calculation_more as TC
from Utils.circuit_cut import list_str_to_int, single_to_all_line_sequence, letter_to_number

from Utils.generate_partitions import generate_line, random_line
from Utils.read_qasm import converter_circ_from_qasm, count_num_of_qubit
from Utils.readreals import real_to_cnot_list, qasm_to_cnot_list, nct_to_ncv

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

# 第二种方法计算全局门数量
'''判断一个门是否为全局门'''
def judge_is_global_gate(gate,cut_point):
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

# 第二种方法计算全局门数量
def count_gg_num(gate_list,cut_point):
    count = 0
    for i in range(len(gate_list)):
        if judge_is_global_gate(gate_list[i],cut_point)==1:
            count+=1
    return count

'''根据线序计算全局门数'''
def count_gg_num_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence):
    new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence)
    gg_num = len(is_global_gate(new_gate_list,cut_list)[1])
    return gg_num

'''根据线序计算全局门数'''
def k_count_gg_num_by_line_sequence(initial_line_sequence,gate_list,line_sequence):
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence)
    gg_num = len(is_global_gate(new_gate_list,cut_list)[1])
    return gg_num

'''根据线序改变gate_list(适用于两个分区)'''
def change_gate_list_by_line_sequence(initial_line_sequence,gate_list,single_line_sequence):
    qubit = max([max(row) for row in gate_list]) + 1
    new_gate_list = list_str_to_int(gate_list)
    all_line_sequence = single_to_all_line_sequence(single_line_sequence, qubit)  # ABCDEG FHIJKL
    # print(all_line_sequence)
    for j in range(len(all_line_sequence)):  # j: 0-11
        if all_line_sequence[j] == initial_line_sequence[j]:
            continue
        if all_line_sequence[j] != initial_line_sequence[j]:
            for k in range(len(new_gate_list)):  # k: 门数
                for p in range(len(new_gate_list[k])):  # l:0-3
                    if new_gate_list[k][p] == letter_to_number(all_line_sequence[j]):
                        new_gate_list[k][p] = str(j)
    new_gate_list = list_str_to_int(new_gate_list)
    return new_gate_list


'''根据线序改变gate_list(适用于K个分区)'''
def k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,all_line_sequence):
    new_gate_list = list_str_to_int(gate_list)
    for j in range(len(all_line_sequence)):  # j: 0-11
        if all_line_sequence[j] == initial_line_sequence[j]:
            continue
        if all_line_sequence[j] != initial_line_sequence[j]:
            for k in range(len(new_gate_list)):  # k: 门数
                for p in range(len(new_gate_list[k])):  # l:0-3
                    if new_gate_list[k][p] == letter_to_number(all_line_sequence[j]):
                        new_gate_list[k][p] = str(j)
    new_gate_list = list_str_to_int(new_gate_list)
    return new_gate_list


'''2分区下计算最小全局门数'''
def count_min_global_gate_num(gate_list,cut_list):
    qubit = max([max(row) for row in gate_list]) + 1
    str1 = ''
    # 最多支持 62量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z']
    for i in range(int(qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    initial_line_sequence = str1
    # 排列组合 将str1按qibt/2 一分为，共有C(qbit/2,qbit)种情况
    line_sequence_combination = []  # 线序排列集合['ABCDEF', 'ABCDEG', 'ABCDEH', 'ABCDEI', 'ABCDEJ', 'ABCDEK', 'ABCDEL', 'ABCDFG', 'ABCDFH', 'ABCDFI', 'ABCDFJ'......]
    min_gg_bum = len(is_global_gate(gate_list,cut_list)[1])
    min_line_sequence = initial_line_sequence[0:cut_list[0]]
    for i in itertools.combinations(str1, cut_list[0]):
        #  print(''.join(i), end=" ")
        # line_sequence_combination.append(''.join(i))
        line_sequence = (''.join(i))
        # all_line_sequence = single_to_all_line_sequence(line_sequence, qubit)  # ABCDEG FHIJKL
        gg_num = count_gg_num_by_line_sequence(initial_line_sequence,gate_list,line_sequence)
        if gg_num < min_gg_bum:
            min_gg_bum = gg_num
            min_line_sequence = line_sequence
    min_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,min_line_sequence)
    return min_gg_bum,min_gate_list


'''K分区下计算最小全局门数'''
def k_count_min_global_gate_num(gate_list,cut_list,line_sequence_list):
    initial_line_sequence = ''
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    qubit = max([max(row) for row in gate_list]) + 1
    for i in range(int(qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        initial_line_sequence = initial_line_sequence + str2[i]
    min_gg_bum = len(is_global_gate(gate_list, cut_list)[1])
    min_line_sequence = initial_line_sequence
    for i in range(len(line_sequence_list)):
        line_sequence = line_sequence_list[i]
        gg_num = count_gg_num_by_line_sequence(initial_line_sequence, gate_list, line_sequence)
        if gg_num < min_gg_bum:
            min_gg_bum = gg_num
            min_line_sequence = line_sequence
    min_gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list, min_line_sequence)
    return min_gg_bum, min_gate_list


'''初始线路'''
def Initial_line_sequence(gate_list):
    initial_line_sequence = ''
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9','$','#']
    qubit = max([max(row) for row in gate_list]) + 1
    for i in range(int(qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        initial_line_sequence = initial_line_sequence + str2[i]
    return initial_line_sequence


'''K分区下计算最小传输代价'''
def k_count_min_st_num2(gate_list,cut_list,line_sequence_list):
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    STO = []
    lenl = len(line_sequence_list)
    print(lenl)
    minio = len(gate_list)*2
    mini = len(gate_list)*2
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end='')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])
        # print(new_gate_list,end=' ')
        # st = TC.direct_calculation_of_tc(new_gate_list, cut_list)
        st,sto = TC.direct_calculation_of_tc_queue(new_gate_list, cut_list)
        ST.append(st)
        STO.append(sto)
        print('当前最小：', minio,' ',mini )
        if i%100 == 0:
            minio = min(STO)
            mini = min(ST)
    return ST

'''K分区下计算最小传输代价'''
def k_count_min_st_num(gate_list,cut_list,line_sequence):

    initial_line_sequence = Initial_line_sequence(gate_list)
    new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence, gate_list, line_sequence)
    # st = TC.direct_calculation_of_tc_distributed(new_gate_list, cut_list,1)
    st,sto = TC.direct_calculation_of_tc_queue(new_gate_list, cut_list)
    return st,sto


def generate_partitions(remaining_letters, partition_sizes, current_partition, all_partitions0,all_partitions1):
    if len(partition_sizes) == 0:
        if len(remaining_letters) == 0:
            line = "".join(["".join(row) for row in current_partition.copy()])
            print(line,end=' ')
            st,sto = k_count_min_st_num(gate_list,cut_list,line)
            all_partitions0.append(st)
            all_partitions1.append(sto)
            print('当前最低传输代价：',end=' ')
            print(min(all_partitions0),min(all_partitions1))
        return

    current_size = partition_sizes[0]
    for comb in itertools.combinations(remaining_letters, current_size):
        next_partition = current_partition.copy()
        next_partition.append(list(comb))
        next_remaining_letters = [letter for letter in remaining_letters if letter not in comb]

        generate_partitions(next_remaining_letters, partition_sizes[1:], next_partition, all_partitions0,all_partitions1)


def generate_line2(partition_sizes):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']

    letters = str2[0:sum(partition_sizes)]

    all_partitions0 = []
    all_partitions1 = []
    generate_partitions(letters, partition_sizes, [], all_partitions0,all_partitions1)

    return all_partitions0

def random_line_by_qubit(random_num,circuit_qubit):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    letters = str2[0:circuit_qubit]
    all_partitions = []
    for i in range(random_num):
        random.shuffle(letters)
        all_partitions.append(''.join(letters))
    return all_partitions




from collections import defaultdict
from typing import List, Tuple


def group_nodes_by_partition(partition_labels: List[int], node_labels: List[str]) -> Tuple[List[int], str]:
    """
    将节点按分区编号分组，返回每个分区的大小和按分区拼接的节点字符串。

    :param partition_labels: 每个节点所属的分区编号
    :param node_labels: 节点名称列表
    :return: (分区大小列表, 拼接后的节点字符串)
    """
    if len(partition_labels) != len(node_labels):
        raise ValueError("partition_labels 和 node_labels 长度必须一致。")

    partition_dict = defaultdict(list)
    for node, partition in zip(node_labels, partition_labels):
        partition_dict[partition].append(node)

    # 按照分区编号排序
    sorted_partitions = sorted(partition_dict.items())

    # 构造结果
    partition_sizes = [len(nodes) for _, nodes in sorted_partitions]
    concatenated_nodes = ''.join([node for _, nodes in sorted_partitions for node in nodes])

    return partition_sizes, concatenated_nodes


# 主函数
if __name__ == '__main__':
    start_time = time.time()
    # 读取qasm
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/reals/new_real/6sym.qasm'
    # gate_list = qasm_to_cnot_list(input_filename)
    # print(gate_list)
    # 读取real
    input_filename = '../reals/new_real/8bitadder.real'
    gate_list = real_to_cnot_list(input_filename)
    # gate_list = nct_to_ncv(input_filename)

    print(gate_list)
    print(len(gate_list))
