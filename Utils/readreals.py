'''
功能：将real文件读取后分解成CNOT门以及单门，并提取CNOT门作为list存储
作者：陈新宇
时间：2023.7.9
'''
from pytket._tket.circuit import Circuit
from pytket._tket.passes import CliffordSimp, KAKDecomposition, RemoveRedundancies
from pytket.qasm import circuit_from_qasm_str, circuit_to_qasm_str
from qiskit import QuantumCircuit, transpile
import re

from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Unroller

from Utils.circuit_cut import list_str_to_int
from Utils.read_qasm import converter_circ_from_qasm


# 前瞻深度
Depth = 1
IS_circuit_opt = 0

'''读取real文件'''
def real_to_list(input_filename):
    real_file = open(input_filename, 'r',encoding='utf-8')

    lines = real_file.readlines()  # 读取所有

    qubit_number = []  # 量子位编号
    for i in range(len(lines)):
        if 'variables' in lines[i]:
            # print(lines[i])
            qubit_number = lines[i][11:-1]
        if 'begin' in lines[i]:  # 获取begin的那一行
            begin = i
    qubit_number_list = qubit_number.split(' ')
    print(qubit_number_list)
    gate_list = lines[begin + 1:len(lines)-1]
    print(gate_list)
    gate_list_all = []
    for i in range(len(gate_list)):
        gg = 0
        gate = []
        fredkin_list = []
        ##去掉\n
        gate_list[i] = gate_list[i].strip()
        e = gate_list[i].split(" ")
        print(e)
        if e[0] == '#':
            continue
        for ind in range(len(e)):
            if ind == 0:
                continue
            if e[0] == 'f3':
                fredkin_gg = qubit_number_list.index(e[ind])
                fredkin_list.append(fredkin_gg)
            else:
                gg = qubit_number_list.index(e[ind])  # 在 qubit_number_list中查找元素位置
                gate.append(gg)

        # fredkin分解 fredkin_list长度为3
        if len(fredkin_list) != 0:
            print(fredkin_list)
            gate_list_all.append([fredkin_list[2], fredkin_list[1]])
            gate_list_all.append([fredkin_list[0], fredkin_list[1], fredkin_list[2]])
            gate_list_all.append([fredkin_list[2], fredkin_list[1]])
        else:
            gate_list_all.append(gate)
    print(gate_list_all)
    return gate_list_all


'''读取tf文件'''
def tf_to_list(input_filename):
    real_file = open(input_filename, 'r',encoding='utf-8')
    lines = real_file.readlines()  # 读取所有
    qubit_number = []  # 量子位编号
    for i in range(len(lines)):
        if '.v' in lines[i]:
            # print(lines[i])
            qubit_number = lines[i][3:-1]
        if 'BEGIN' in lines[i]:  # 获取begin的那一行
            begin = i
    qubit_number_list = qubit_number.split(' ')
    print(qubit_number_list)
    gate_list = lines[begin + 1:len(lines)-1]
    print(gate_list)
    gate_list_all = []
    for i in range(len(gate_list)):
        gate = []
        ##去掉\n
        gate_list[i] = gate_list[i].strip()
        e = gate_list[i].split(" ")
        print(e)
        if e[0] == '#':
            continue
        for ind in range(len(e)):
            if ind == 0:
                continue
            else:
                gg = qubit_number_list.index(e[ind])  # 在 qubit_number_list中查找元素位置
                gate.append(gg)
        gate_list_all.append(gate)
    print(gate_list_all)
    print(len(gate_list_all))
    return gate_list_all


'''gate_list转circuit'''
def gate_list_to_circuit(gate_list):
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    circuit = QuantumCircuit(circuit_qubit, circuit_qubit)
    for i in range(len(gate_list)):
        if len(gate_list[i]) == 1:
            circuit.x(gate_list[i][0])
        if len(gate_list[i]) == 2:
            circuit.cx(gate_list[i][0], gate_list[i][1])
        if len(gate_list[i]) >= 3:
            control_bits = gate_list[i][0:-1]
            circuit.mcx(control_bits, gate_list[i][-1])
    # print(circuit)
    # decomposed_circ = circuit.decompose()
    # decomposed_circ = transpile(circuit, basis_gates=['u1', 'u2', 'u3', 'cx'], optimization_level=1)
    # CliffordSimp().apply(circuit)
    # print(circuit)
    # decomposed_circ = circuit
    # pass_ = Unroller(['u1', 'u2', 'u3', 'cx'])
    # pm = PassManager(pass_)
    # decomposed_circ = pm.run(circuit)
    # print(decomposed_circ)
    # print(555)
    # ops_list = [(op.name, qargs) for op, qargs, cargs in decomposed_circ.data]
    # # Print the list of operations
    # for j in range(len(ops_list)):
    #     print(ops_list[j])
    return circuit

'''gate_list转circuit 一个个门转'''
def gate_list_to_circuit_tket(gate_list):
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    circ = QuantumCircuit(circuit_qubit)
    gate_count = 0
    for i in range(len(gate_list)):
        print(gate_list[i])
        if len(gate_list[i]) == 1:
            circ.x(gate_list[i][0])
            gate_count += 1
        if len(gate_list[i]) == 2:
            circ.cx(gate_list[i][0], gate_list[i][1])
            gate_count += 2
        if len(gate_list[i]) >= 3:
            # 对单个MCT门分解成ggate_list
            control_bits = gate_list[i][0:-1]
            target_bits = gate_list[i][-1]
            qcircuit = QuantumCircuit(circuit_qubit)
            qcircuit.mct(control_bits, target_bits)
            decomposed_circ = transpile(qcircuit, basis_gates=['u1', 'u2', 'u3', 'cx'], optimization_level=3)
            gate_count += sum(decomposed_circ.count_ops().values())
            ggate_list = qasm_list_to_gate_lsit(qc_to_qasm_list(decomposed_circ))
            for j in range(len(ggate_list)):
                if len(ggate_list[j]) == 2:
                    circ.cx(int(ggate_list[j][0]), int(ggate_list[j][1]))

    print('量子门数：',end=' ')
    print(gate_count)
    return circ


def split_list(input_list, d):
    return [input_list[i:i+d] for i in range(0, len(input_list), d)]

'''gate_list转circuit '''
def gate_list_to_circuit_tket_depth(gate_list,depth):
    circuit_qubit = max([max(row) for row in gate_list]) + 10
    circuit = QuantumCircuit(circuit_qubit)
    sub_gate_list = split_list(gate_list, depth)
    print(sub_gate_list)
    for i in range(len(sub_gate_list)):
        print(sub_gate_list[i])
        sub_circ = gate_list_to_circuit(sub_gate_list[i])
        decomposed_circ = transpile(sub_circ, basis_gates=['h', 'rx', 'ry', 'rz', 'x', 'y', 'z', 's', 't', 'tdg', 'sdg', 'cz','cx'], optimization_level=3)
        # print(decomposed_circ)
        # decomposed_circ = transpile(sub_circ,
        #                                 basis_gates=['h','u1','u2','u3','cz'],
        #                                 optimization_level=3)
        # decomposed_circ = transpile(sub_circ, basis_gates=['h', 'u1', 'u2','cu1', 'cx'], optimization_level=3)
        # print(decomposed_circ.qasm())
        if IS_circuit_opt == 1:
            circ = circuit_from_qasm_str(decomposed_circ.qasm())
            CliffordSimp().apply(circ)
            RemoveRedundancies().apply(circ)
            decomposition_pass = KAKDecomposition()
            decomposition_pass.apply(circ)
            circuit_opt_qasm = circuit_to_qasm_str(circ)
            decomposed_circ = QuantumCircuit.from_qasm_str(circuit_opt_qasm)


        ggate_list = qasm_list_to_gate_lsit(qc_to_qasm_list(decomposed_circ))
        for j in range(len(ggate_list)):
            if len(ggate_list[j]) == 2:
                circuit.cx(int(ggate_list[j][0]), int(ggate_list[j][1]))
    return circuit



def qc_to_qasm_list(qc):
    # qc 转化为 qsam
    qc_qasm = qc.qasm()  # str类型
    # print(qc_qasm)
    # print(qc_qasm.count(';'))
    # print(qc_qasm.find(';'))
    lst = []
    # ；的索引号，用于换行
    for pos, char in enumerate(qc_qasm):
        if (char == ';'):
            lst.append(pos)
    # print(lst)
    qc_qasm_list = []
    qc_qasm_list.append(qc_qasm[0:lst[0] + 1])
    # print(qc_qasm_list)
    for i in range(len(lst) - 1):
        qc_qasm_list.append(qc_qasm[lst[i] + 2:lst[i + 1] + 1])
    # print(qc_qasm_list)
    return qc_qasm_list

''' 获取数字 '''
def get_data(str):
    pattern = re.compile(r'\d+')  # 查找数字
    result = pattern.findall(str)
    return result

''' 将list格式的qasm转化为datefrome格式 '''
def qasm_list_to_gate_lsit(qc_qasm_list):
    gate_list = []
    for i in range(4, len(qc_qasm_list)):  # 一行行遍历
        line = qc_qasm_list[i]
        # print(line)
        if line[0:2] == 'CX' or line[0:2] == 'cx':
            '''获取CNOT'''
            cnot = get_data(line)
            cnot_control = cnot[0]
            cnot_target = cnot[1]
            gate_list.append([cnot_control,cnot_target])
        if line[0:2] == 'CZ' or line[0:2] == 'cz':
            '''获取CNOT'''
            cnot = get_data(line)
            cnot_control = cnot[0]
            cnot_target = cnot[1]
            gate_list.append([cnot_control,cnot_target])
        if line[0:3] == 'cu1':
            '''获取Cu'''
            cnot = get_data(line)
            cnot_control = cnot[0]
            cnot_target = cnot[1]
            gate_list.append([cnot_control, cnot_target])
    return gate_list




def real_to_cnot_list(input_filename):
    # fredkin
    # Circuit = QuantumCircuit(6,6)
    # Circuit.fredkin(0,2,3)
    # Circuit.fredkin(1, 2, 4)
    # Circuit.fredkin(1, 3, 5)
    # decomposed_circ = Circuit.decompose()
    # 其他
    gate_list = real_to_list(input_filename)
    # circuit = gate_list_to_circuit_tket_depth(gate_list,Depth)
    circuit = gate_list_to_circuit_tket_depth(gate_list,len(gate_list))
    print(circuit.qasm())
    print('量子门数：',end=' ')
    print(circuit.count_ops())
    # gate_count = circuit.count_ops()
    # for gate, count in gate_count.items():
    #     print(f"{gate}: {count}")
    print(666)
    qc_qasm_list = qc_to_qasm_list(circuit)
    print(777)
    gate_list = qasm_list_to_gate_lsit(qc_qasm_list)
    print(888)
    return list_str_to_int(gate_list)

def real_to_qasm(input_filename):
    gate_list = real_to_list(input_filename)
    circuit = gate_list_to_circuit(gate_list)
    qasm = circuit.qasm()
    return qasm

def nct_to_ncv(input_filename):
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    ncv_list = []
    for gate in gate_list:
        if len(gate) == 2:
            ncv_list.append(gate)
        if len(gate) == 3:
            ncv_list.append([gate[1],gate[2]])
            ncv_list.append([gate[0], gate[1]])
            ncv_list.append([gate[1], gate[2]])
            ncv_list.append([gate[0], gate[1]])
            ncv_list.append([gate[0], gate[2]])
    return ncv_list

def nct_to_ncv_2(gate_list):
    ncv_list = []
    for gate in gate_list:
        if len(gate) == 2:
            ncv_list.append(gate)
        if len(gate) == 3:
            ncv_list.append([gate[1],gate[2]])
            ncv_list.append([gate[0], gate[1]])
            ncv_list.append([gate[1], gate[2]])
            ncv_list.append([gate[0], gate[1]])
            ncv_list.append([gate[0], gate[2]])
    return ncv_list


'''全过程'''
def qasm_to_cnot_list(input_filename):
    # fredkin
    # Circuit = QuantumCircuit(6,6)
    # Circuit.fredkin(0,2,3)
    # Circuit.fredkin(1, 2, 4)
    # Circuit.fredkin(1, 3, 5)
    # decomposed_circ = Circuit.decompose()
    # 其他
    gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    print(gate_list)
    # circuit = gate_list_to_circuit_tket(gate_list)
    circuit = gate_list_to_circuit_tket_depth(gate_list, len(gate_list))
    gate_count = circuit.count_ops()
    print(gate_count)
    qc_qasm_list = qc_to_qasm_list(circuit)
    gate_list = qasm_list_to_gate_lsit(qc_qasm_list)
    return list_str_to_int(gate_list)


# 统计空量子位
def count_empty_qubit(gate):
    # gate = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 20, 17]
    empty_qubit = max(gate)-min(gate)-len(gate)+1
    return empty_qubit

def compose_rule2(gate):
    gate_list = []
    n = gate[-1]-gate[0]+1
    control_qubits_num = gate[0:-1]
    gate1_control_qubits = int(len(gate) / 2) + 1
    gate2_control_qubits = n-1-gate1_control_qubits
    print(gate1_control_qubits)
    print(gate2_control_qubits)
    gate1 = gate[0:gate1_control_qubits]
    gate1.append()
    return gate_list

# 递归分解
def decompose_mct(gate,gate_list=None):
    if gate_list is None:
        gate_list = []
    if len(gate) <= 2:
        gate_list.append(gate)
    if len(gate) > 2:
        gate_list.append([gate[-2], gate[-1]])
        decompose_mct(gate[0:-1],gate_list)
        gate_list.append([gate[-2], gate[-1]])
        decompose_mct(gate[0:-1],gate_list)
        v_gate = gate[0:-2]
        v_gate.append(gate[-1])
        decompose_mct(v_gate,gate_list)

    return gate_list

def resl_to_decompose_ncv(input_filename):
    gate_list = real_to_list(input_filename)
    ncv_list = []
    for gate in gate_list:
        ncv_list += decompose_mct(gate)
    return ncv_list


if __name__ == '__main__':
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/reals/cycle17_3_112.real'
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/reals/urf6_160.real'
    # gate_list = resl_to_decompose_ncv(input_filename)
    # gate_list = real_to_cnot_list(input_filename)

    # gate_listg = real_to_list(input_filename)
    # gate_listg.reverse()
    # print(gate_listg)

    # for k in range(len(gate_listg)):
    #     gate_list = [gate_listg[k]]
    #     print(gate_list)
    #
    #     circuit = gate_list_to_circuit(gate_list)
    #     print(circuit)
    #     decompose_circuit = transpile(circuit, basis_gates=['h', 'u1', 'cu1', 'cx'], optimization_level=3)
    #     count = 0
    #     ggate_list = qasm_list_to_gate_lsit(qc_to_qasm_list(decompose_circuit))
    #     for j in range(len(ggate_list)):
    #         if len(ggate_list[j]) == 2:
    #             count += 1
    #     print(count)


    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/reals/new_real/cycle17_3-corrected.qc'
    # gate_list = tf_to_list(input_filename)



    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/8bitadder2.qasm'
    # # gate_list = nct_to_ncv(input_filename)
    # gate_list = qasm_to_cnot_list(input_filename)
    # print('cnot门数：',end=' ')
    # print(len(gate_list))


    # circ = gate_list_to_circuit_tket_depth(gate_list,1)
    # print(circ.qasm())

    input_filename = '../reals/alu_318.real'
    gate_list = real_to_cnot_list(input_filename)
    # print(gate_list)
    print(len(gate_list))
    # 另存为txt
    # 将列表写入文本文件
    file_path = "../Txt/alu_318.txt"  # 文件路径

    with open(file_path, "w") as file:
        for sublist in gate_list:
            file.write(" ".join(map(str, sublist)) + "\n")




    # # # 另存为qasm
    # # gate_list = real_to_list(input_filename)
    # circuit = gate_list_to_circuit(gate_list)
    # qasm = circuit.qasm()
    # # print(qasm)


    # qasm_file_path = "../synthesis/15_and_16_qubits_test/16qubit_circuit/hwb_12.qasm"
    #
    # with open(qasm_file_path, "w") as qasm_file:
    #     qasm_file.write(qasm)
