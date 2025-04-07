'''
功能： 代码主要面向于线路划分为三部分,四部分
时间： 2023.03.08
版本： 1.0
作者： 陈新宇
'''

import re
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import Utils.write_txt as txt
import itertools


from Utils.draw_circuit import new_circuit
from Utils.read_qasm import count_num_of_qubit, remove_single_qubit_gate
from Utils.transmission_cost_calculation_2 import taboo_search
from Utils.transmission_cost_calculation_more import direct_calculation_of_tc, is_global_gate, statistics_gate_labels, \
    transfer_qubit_list_by_gate_list, count_transfer_queue



'''读取qasm文件并进行存储'''
def converter_circ_from_qasm(input_file_name):
    gate_list = []
    qbit = 0  # 量子位
    qasm_file = open(input_file_name, 'r')
    iter_f = iter(qasm_file)
    reserve_line = 0
    num_line = 0
    for line in iter_f:  # 遍历文件，一行行遍历，读取文本
        num_line += 1
        if num_line <= reserve_line:
            continue
        else:
            if 'qreg' in line:
                qbit = get_data(line)[0]
            if line[0:1] == 'x' or line[0:1] == 'X':
                '''获取X门'''
                x = get_data(line)
                x_target = x[0]
                listSingle = [x_target]
                gate_list.append(listSingle)
            if line[0:2] == 'CX' or line[0:2] == 'cx':
                '''获取CNOT'''
                cnot = get_data(line)
                cnot_control = cnot[0]
                cnot_target = cnot[1]
                listSingle = [cnot_control, cnot_target]
                gate_list.append(listSingle)
            if line[0:2] == 'CP' or line[0:2] == 'cp':
                cp = get_data(line)
                cp_one = cp[1]
                cp_two = cp[2]
                listSingle = [cp_one,cp_two]
                gate_list.append(listSingle)
            if line[0:4] == 'SWAP' or line[0:4] == 'swap':
                swap = get_data(line)
                swap_one = swap[0]
                swap_two = swap[1]
                cnot_one = [swap_one,swap_two]
                cnot_two = [swap_two,swap_one]
                gate_list.append(cnot_one)
                gate_list.append(cnot_two)
                gate_list.append(cnot_one)
            if line[0:3] == 'CCX' or line[0:3] == 'ccx':
                '''获取toffoli'''
                toffoli = get_data(line)
                toffoli_control1 = toffoli[0]
                toffoli_control2 = toffoli[1]
                toffoli_target = toffoli[2]
                listSingle = [toffoli_control1, toffoli_control2, toffoli_target]
                gate_list.append(listSingle)
    return gate_list, qbit

def get_data(str):
    pattern = re.compile("[\d]+")
    result = re.findall(pattern, str)
    return result

'''
将gate_list全部转换为int
'''
def list_str_to_int(gate_list):
    new_gate_list = []
    for i in range(len(gate_list)):
        son_new_gate_list = list(map(int, gate_list[i]))
        new_gate_list.append(son_new_gate_list)
    return new_gate_list



'''
根据前半部分线序推出全部线序  ABCDFG =>ABCDFGEHIJKL
'''
def single_to_all_line_sequence(line_sequence_combination,qbit):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    all_line_sequence = list(line_sequence_combination)
    for i in range(int(qbit)):
        if str2[i] not in all_line_sequence:
            all_line_sequence.append(str2[i])
    new_all_line_sequence = "".join(all_line_sequence)
    return new_all_line_sequence

'''
根据前半部分线序推出后半线序  ABCDFG =>EHIJKL
'''
def single_to_leave_line_sequence(line_sequence_combination,qbit):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    all_line_sequence = list(line_sequence_combination)
    leave_line_sequence = []
    for i in range(int(qbit)):
        if str2[i] not in all_line_sequence:
            leave_line_sequence.append(str2[i])
    new_leave_line_sequence = "".join(leave_line_sequence)
    return new_leave_line_sequence

'''
二分法下输出第一次划分所有种线序组合 C(cut_point,qbit)
'''
def line_sequence_change_combination_2(qbit,cut_list):
    str1 = ''
    cut_point = cut_list[0]
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(qbit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    # 排列组合 将str1按qibt/2 一分为，共有C(qbit/2,qbit)种情况
    line_sequence_combination = []  # 线序排列集合['ABCDEF', 'ABCDEG', 'ABCDEH', 'ABCDEI', 'ABCDEJ', 'ABCDEK', 'ABCDEL', 'ABCDFG', 'ABCDFH', 'ABCDFI', 'ABCDFJ'......]
    for i in itertools.combinations(str1, cut_point):
        #  print(''.join(i), end=" ")
        # line_sequence_combination.append(''.join(i))
        print(single_to_all_line_sequence(''.join(i),qbit))
        line_sequence_combination.append(single_to_all_line_sequence(''.join(i),qbit))
    return line_sequence_combination


'''三分法下线序'''
def line_sequence_change_combination_3(qbit,cut_list):
    str1 = ''
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(qbit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    # 排列组合
    line_sequence_combination = []  # 线序排列集合['ABCDEF', 'ABCDEG', 'ABCDEH', 'ABCDEI', 'ABCDEJ', 'ABCDEK', 'ABCDEL', 'ABCDFG', 'ABCDFH', 'ABCDFI', 'ABCDFJ'......]
    for i in itertools.combinations(str1, cut_list[0]):
        line_sequence_1 = ''.join(i)
        line_sequence_23 = single_to_leave_line_sequence(line_sequence_1,qbit)
        # 对23部分线序重新组合
        for j in itertools.combinations(line_sequence_23, cut_list[1]):
            line_sequence_2 = ''.join(j)
            line_sequence_12 = line_sequence_1+line_sequence_2
            # print(line_sequence_2,end=' ')
            line_sequence_3 = single_to_leave_line_sequence(line_sequence_12,qbit)
            line_sequence = line_sequence_1+line_sequence_2+line_sequence_3
            line_sequence_combination.append(line_sequence)
    return line_sequence_combination


'''四分法下线序'''
def line_sequence_change_combination_4(qbit,cut_list):
    str1 = ''
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(qbit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    # 排列组合
    line_sequence_combination = []  # 线序排列集合['ABCDEF', 'ABCDEG', 'ABCDEH', 'ABCDEI', 'ABCDEJ', 'ABCDEK', 'ABCDEL', 'ABCDFG', 'ABCDFH', 'ABCDFI', 'ABCDFJ'......]
    for i in itertools.combinations(str1, cut_list[0]):
        line_sequence_1 = ''.join(i)
        line_sequence_234 = single_to_leave_line_sequence(line_sequence_1,qbit)
        # 对23部分线序重新组合
        for j in itertools.combinations(line_sequence_234, cut_list[1]):
            line_sequence_2 = ''.join(j)
            line_sequence_12 = line_sequence_1+line_sequence_2
            line_sequence_34 = single_to_leave_line_sequence(line_sequence_12, qbit)
            for k in itertools.combinations(line_sequence_34, cut_list[2]):
                line_sequence_3 = ''.join(k)
                line_sequence_123 = line_sequence_12 + line_sequence_3
                line_sequence_4 = single_to_leave_line_sequence(line_sequence_123, qbit)
                line_sequence = line_sequence_123+line_sequence_4
                line_sequence_combination.append(line_sequence)
    return line_sequence_combination

'''生成初始线序'''
def new_initial_line_sequence(circuit_qubit):
    initial_line_sequence = ''
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(circuit_qubit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        initial_line_sequence = initial_line_sequence + str2[i]
    return initial_line_sequence


'''根据线序改变gate_list'''
def change_gate_list_by_line_sequence(initial_line_sequence,gate_list,all_line_sequence):
    new_gate_list = list_str_to_int(gate_list)
    circuit_qbit =  max([max(row) for row in new_gate_list]) + 1
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

'''
把字母变成对于的数字
'''
def letter_to_number(letter):
    number = ord(letter) - 65
    return number


'''线路分割成两部分,计算最低tc'''
def circuit_cut_2(gate_list,cut_list):
    str1 = ''
    qbit = max([max(row) for row in gate_list]) + 1
    cut_point = cut_list[0]
    initial_line_sequence = new_initial_line_sequence(qbit)
    # 最多支持 26量子位
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z']
    for i in range(int(qbit)):
        # str1：ABCDEFGHIJKL.... 共qbit个
        str1 = str1 + str2[i]
    # 排列组合 将str1按qibt/2 一分为，共有C(qbit/2,qbit)种情况
    best_line_sequence_combination = []
    min_tc = direct_calculation_of_tc(gate_list,cut_list)
    for i in itertools.combinations(str1, cut_point):
        line_sequence = single_to_all_line_sequence(''.join(i), qbit)
        print(line_sequence,end=' ')
        new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence)
        new_tc = direct_calculation_of_tc(new_gate_list,cut_list)
        print(new_tc)
        if new_tc == min_tc:
            best_line_sequence_combination.append(line_sequence)
        if new_tc < min_tc:
            min_tc = new_tc
            best_line_sequence_combination.clear()
            best_line_sequence_combination.append(line_sequence)

    return min_tc,best_line_sequence_combination




if __name__ == '__main__':
    # 读取线路
    input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/mini-alu_167.qasm'
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename)[0])
    gate_list = remove_single_qubit_gate(gate_list)
    gate_list = [[4, 1], [2, 4], [0, 2], [4, 0], [4, 2], [0, 2], [4, 0], [2, 4], [5, 0], [4, 5], [0, 4], [0, 5], [4, 5], [0, 4], [5, 0], [2, 3], [5, 2], [3, 5], [3, 2], [5, 2], [3, 5], [2, 3], [5, 0], [4, 5], [0, 4], [0, 5], [4, 5], [0, 4], [5, 0], [2, 3], [5, 2], [3, 5], [3, 2], [5, 2], [3, 5], [2, 3], [0, 1], [4, 0], [1, 4], [1, 0], [4, 0], [1, 4], [0, 1], [2, 3], [5, 2], [3, 5], [3, 2], [5, 2], [3, 5], [2, 3], [3, 5], [4, 3], [5, 4], [5, 3], [4, 3], [5, 4], [3, 5], [0, 1], [3, 0], [1, 3], [1, 0], [3, 0], [1, 3], [0, 1], [3, 5], [4, 3], [5, 4], [5, 3], [4, 3], [5, 4], [3, 5], [0, 1], [3, 0], [1, 3], [1, 0], [3, 0], [1, 3], [0, 1], [2, 3], [5, 2], [3, 5], [3, 2], [5, 2], [3, 5], [2, 3], [3, 5], [4, 3], [5, 4], [5, 3], [4, 3], [5, 4], [3, 5], [0, 1], [3, 0], [1, 3], [1, 0], [3, 0], [1, 3], [0, 1], [3, 5], [4, 3], [5, 4], [5, 3], [4, 3], [5, 4], [3, 5], [0, 1], [3, 0], [1, 3], [1, 0], [3, 0], [1, 3], [0, 1]]


    print('qasm读取线路:', end=' ')
    print(gate_list)
    qubit_num = count_num_of_qubit(gate_list)
    print('量子位数：' + str(qubit_num))

    cut_list = [3,3]

    # print(line_sequence_change_combination_3(4,[1,1,2]))
    # print(len(line_sequence_change_combination_4(4,[1,1,1,1])))

    min_tc,best_line_sequence_combination = circuit_cut_2(gate_list, cut_list)
    print(min_tc)
    print(best_line_sequence_combination)

    initial_line_sequence = new_initial_line_sequence(qubit_num)
    new_gate_list = change_gate_list_by_line_sequence(initial_line_sequence, gate_list,best_line_sequence_combination[0])
    print(new_gate_list)


    gate_list = new_gate_list
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
    print('全局门：' + str(global_gate_num_list))
    print('初始传输队列：' + str(initial_transfer_queue))
    print('初始传输代价：' + str(len(initial_transfer_queue) * 2))
    print('######################################################################################################')









