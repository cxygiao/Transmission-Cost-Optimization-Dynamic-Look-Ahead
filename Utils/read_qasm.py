
import re

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
            if line[0:2] == 'CX' or line[0:2] == 'cx' or line[0:2] == 'cz':
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
    gate_list = remove_single_qubit_gate(gate_list)
    return gate_list

def get_data(str):
    pattern = re.compile("[\d]+")
    result = re.findall(pattern, str)
    return result

'''获取量子位数'''
def count_num_of_qubit(gate_list):
    # 线路的量子位数
    num_of_qubit = max([max(row) for row in gate_list]) + 1
    return num_of_qubit

'''gate_list去除单门'''
def remove_single_qubit_gate(gate_list):
    i = 0
    while i<len(gate_list):
        if len(gate_list[i]) == 1 :
            gate_list.pop(i)
            i -= 1
        i += 1
    return gate_list

'''
将gate_list全部转换为int
'''
def list_str_to_int(gate_list):
    new_gate_list = []
    for i in range(len(gate_list)):
        son_new_gate_list = list(map(int, gate_list[i]))
        new_gate_list.append(son_new_gate_list)
    return new_gate_list