from qiskit import QuantumRegister, ClassicalRegister, AncillaRegister, QuantumCircuit, transpiler, transpile
from qiskit.circuit.library import MCXGate
import time, os
from PIL import Image
import qiskit.circuit.library.basis_change.qft as qft
from qiskit.circuit import Parameter

# list生成电路图
from Utils.read_qasm import list_str_to_int, converter_circ_from_qasm


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
    circuit = transpile(circuit)
    # 绘制线路
    circuit.draw('mpl', scale=1, filename=img_path, initial_state=True, plot_barriers=False, justify='none', fold=-1)
    # circuit.draw('mpl',scale=1)
    print(circuit)
    return img_path

# list生成电路图
def new_circuit2(gate_list):
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
    return circuit

#circuit生成电路图
def new_qft_circuit(qubit):
    circ = qft.QFT(qubit)
    print(circ)
    circ.draw('mpl', scale=1, filename=img_path, initial_state=True, plot_barriers=False, justify='left', fold=-1)
    return img_path


if __name__ == '__main__':
    # 解决 pillow 加载超大图片报错问题
    Image.MAX_IMAGE_PIXELS = None
    Image.LOAD_TRUNCATED_IMAGES = True

    img_name = time.strftime("%Y-%m-%d %H.%M.%S", time.localtime()) + '.png'
    img_path = os.path.join(os.getcwd(), 'PNG', img_name)
    if not os.path.exists(os.path.join(os.getcwd(), 'PNG')):
        os.mkdir(os.path.join(os.getcwd(), 'PNG'))


    # gate_list = [[0, 1, 9], [2, 9, 10], [1, 2, 9], [3, 10, 11], [3, 9, 10], [2, 3, 9], [4, 10, 11], [4, 9, 10], [3, 4, 9], [6, 10, 11], [6, 9, 10], [4, 6, 9], [7, 10, 11], [7, 9, 10], [6, 7, 9], [8, 10, 11], [8, 9, 10], [7, 8, 9], [5, 10, 11], [5, 9, 10], [10, 11]]

    # gate_list = [[6, 7, 2], [1, 2, 3], [6, 7, 8], [7, 8], [8], [0, 8, 3], [0, 3, 8], [6, 9], [7, 9], [1, 4], [2, 4], [1, 9, 4], [1, 2, 4], [1, 5], [1, 2, 5], [2, 5], [5], [0, 5, 4], [0, 4, 5], [5]]
    # gate_list = [[1, 2, 4], [3, 4, 5], [1, 2, 6], [2, 6], [6], [0, 6, 5], [0, 5, 6], [1, 7], [2, 7], [3, 8], [4, 8], [3, 7, 8], [3, 4, 8], [3, 9], [3, 4, 9], [4, 9], [9], [0, 9, 8], [0, 8, 9], [9]]
    # gate_list = [[31], [31, 30], [31, 29], [31, 28], [31, 27], [31, 26], [31, 9], [31, 8], [31, 31], [31, 30], [31, 29], [31, 28], [31, 27], [31, 26], [31, 25], [31, 24], [31, 23], [31, 22], [31, 21], [31, 20], [31, 19], [31, 18], [31, 17], [31, 16], [31, 7], [31, 6], [31, 5], [31, 4], [31, 3], [31, 2], [31, 1], [31, 0], [30], [30, 29], [30, 28], [30, 27], [30, 26], [30, 9], [30, 8], [30, 31], [30, 30], [30, 29], [30, 28], [30, 27], [30, 26], [30, 25], [30, 24], [30, 23], [30, 22], [30, 21], [30, 20], [30, 19], [30, 18], [30, 17], [30, 16], [30, 7], [30, 6], [30, 5], [30, 4], [30, 3], [30, 2], [30, 1], [30, 0], [29], [29, 28], [29, 27], [29, 26], [29, 9], [29, 8], [29, 31], [29, 30], [29, 29], [29, 28], [29, 27], [29, 26], [29, 25], [29, 24], [29, 23], [29, 22], [29, 21], [29, 20], [29, 19], [29, 18], [29, 17], [29, 16], [29, 7], [29, 6], [29, 5], [29, 4], [29, 3], [29, 2], [29, 1], [29, 0], [28], [28, 27], [28, 26], [28, 9], [28, 8], [28, 31], [28, 30], [28, 29], [28, 28], [28, 27], [28, 26], [28, 25], [28, 24], [28, 23], [28, 22], [28, 21], [28, 20], [28, 19], [28, 18], [28, 17], [28, 16], [28, 7], [28, 6], [28, 5], [28, 4], [28, 3], [28, 2], [28, 1], [28, 0], [27], [27, 26], [27, 9], [27, 8], [27, 31], [27, 30], [27, 29], [27, 28], [27, 27], [27, 26], [27, 25], [27, 24], [27, 23], [27, 22], [27, 21], [27, 20], [27, 19], [27, 18], [27, 17], [27, 16], [27, 7], [27, 6], [27, 5], [27, 4], [27, 3], [27, 2], [27, 1], [27, 0], [26], [26, 9], [26, 8], [26, 31], [26, 30], [26, 29], [26, 28], [26, 27], [26, 26], [26, 25], [26, 24], [26, 23], [26, 22], [26, 21], [26, 20], [26, 19], [26, 18], [26, 17], [26, 16], [26, 7], [26, 6], [26, 5], [26, 4], [26, 3], [26, 2], [26, 1], [26, 0], [9], [9, 8], [9, 31], [9, 30], [9, 29], [9, 28], [9, 27], [9, 26], [9, 25], [9, 24], [9, 23], [9, 22], [9, 21], [9, 20], [9, 19], [9, 18], [9, 17], [9, 16], [9, 7], [9, 6], [9, 5], [9, 4], [9, 3], [9, 2], [9, 1], [9, 0], [8], [8, 31], [8, 30], [8, 29], [8, 28], [8, 27], [8, 26], [8, 25], [8, 24], [8, 23], [8, 22], [8, 21], [8, 20], [8, 19], [8, 18], [8, 17], [8, 16], [8, 7], [8, 6], [8, 5], [8, 4], [8, 3], [8, 2], [8, 1], [8, 0], [31], [31, 30], [31, 29], [31, 28], [31, 27], [31, 26], [31, 25], [31, 24], [31, 23], [31, 22], [31, 21], [31, 20], [31, 19], [31, 18], [31, 17], [31, 16], [31, 7], [31, 6], [31, 5], [31, 4], [31, 3], [31, 2], [31, 1], [31, 0], [30], [30, 29], [30, 28], [30, 27], [30, 26], [30, 25], [30, 24], [30, 23], [30, 22], [30, 21], [30, 20], [30, 19], [30, 18], [30, 17], [30, 16], [30, 7], [30, 6], [30, 5], [30, 4], [30, 3], [30, 2], [30, 1], [30, 0], [29], [29, 28], [29, 27], [29, 26], [29, 25], [29, 24], [29, 23], [29, 22], [29, 21], [29, 20], [29, 19], [29, 18], [29, 17], [29, 16], [29, 7], [29, 6], [29, 5], [29, 4], [29, 3], [29, 2], [29, 1], [29, 0], [28], [28, 27], [28, 26], [28, 25], [28, 24], [28, 23], [28, 22], [28, 21], [28, 20], [28, 19], [28, 18], [28, 17], [28, 16], [28, 7], [28, 6], [28, 5], [28, 4], [28, 3], [28, 2], [28, 1], [28, 0], [27], [27, 26], [27, 25], [27, 24], [27, 23], [27, 22], [27, 21], [27, 20], [27, 19], [27, 18], [27, 17], [27, 16], [27, 7], [27, 6], [27, 5], [27, 4], [27, 3], [27, 2], [27, 1], [27, 0], [26], [26, 25], [26, 24], [26, 23], [26, 22], [26, 21], [26, 20], [26, 19], [26, 18], [26, 17], [26, 16], [26, 7], [26, 6], [26, 5], [26, 4], [26, 3], [26, 2], [26, 1], [26, 0], [25], [25, 24], [25, 23], [25, 22], [25, 21], [25, 20], [25, 19], [25, 18], [25, 17], [25, 16], [25, 7], [25, 6], [25, 5], [25, 4], [25, 3], [25, 2], [25, 1], [25, 0], [24], [24, 23], [24, 22], [24, 21], [24, 20], [24, 19], [24, 18], [24, 17], [24, 16], [24, 7], [24, 6], [24, 5], [24, 4], [24, 3], [24, 2], [24, 1], [24, 0], [23], [23, 22], [23, 21], [23, 20], [23, 19], [23, 18], [23, 17], [23, 16], [23, 7], [23, 6], [23, 5], [23, 4], [23, 3], [23, 2], [23, 1], [23, 0], [22], [22, 21], [22, 20], [22, 19], [22, 18], [22, 17], [22, 16], [22, 7], [22, 6], [22, 5], [22, 4], [22, 3], [22, 2], [22, 1], [22, 0], [21], [21, 20], [21, 19], [21, 18], [21, 17], [21, 16], [21, 7], [21, 6], [21, 5], [21, 4], [21, 3], [21, 2], [21, 1], [21, 0], [20], [20, 19], [20, 18], [20, 17], [20, 16], [20, 7], [20, 6], [20, 5], [20, 4], [20, 3], [20, 2], [20, 1], [20, 0], [19], [19, 18], [19, 17], [19, 16], [19, 7], [19, 6], [19, 5], [19, 4], [19, 3], [19, 2], [19, 1], [19, 0], [18], [18, 17], [18, 16], [18, 7], [18, 6], [18, 5], [18, 4], [18, 3], [18, 2], [18, 1], [18, 0], [17], [17, 16], [17, 7], [17, 6], [17, 5], [17, 4], [17, 3], [17, 2], [17, 1], [17, 0], [16], [16, 7], [16, 6], [16, 5], [16, 4], [16, 3], [16, 2], [16, 1], [16, 0], [7], [7, 6], [7, 5], [7, 4], [7, 3], [7, 2], [7, 1], [7, 0], [6], [6, 5], [6, 4], [6, 3], [6, 2], [6, 1], [6, 0], [5], [5, 4], [5, 3], [5, 2], [5, 1], [5, 0], [4], [4, 3], [4, 2], [4, 1], [4, 0], [3], [3, 2], [3, 1], [3, 0], [2], [2, 1], [2, 0], [1], [1, 0], [0], [0, 31], [0, 31], [0, 31], [1, 30], [1, 30], [1, 30], [2, 29], [2, 29], [2, 29], [3, 28], [3, 28], [3, 28], [4, 27], [4, 27], [4, 27], [5, 26], [5, 26], [5, 26], [6, 9], [6, 9], [6, 9], [7, 8], [7, 8], [7, 8], [16, 31], [16, 31], [16, 31], [17, 30], [17, 30], [17, 30], [18, 29], [18, 29], [18, 29], [19, 28], [19, 28], [19, 28], [20, 27], [20, 27], [20, 27], [21, 26], [21, 26], [21, 26], [22, 25], [22, 25], [22, 25], [23, 24], [23, 24], [23, 24]]
    # gate_list = [[2, 1], [4, 2], [1, 4], [1, 2], [4, 2], [1, 4], [2, 1], [4, 3], [5, 4], [3, 5], [3, 4], [5, 4], [3, 5],
    #              [4, 3], [2, 1], [6, 2], [1, 6], [1, 2], [6, 2], [1, 6], [2, 1], [2, 6], [6, 0], [5, 6], [0, 5], [0, 6],
    #              [5, 6], [0, 5], [6, 0], [5, 0], [6, 5], [0, 6], [0, 5], [6, 5], [0, 6], [5, 0], [1, 7], [2, 7], [3, 8],
    #              [4, 8], [7, 3], [8, 7], [3, 8], [3, 7], [8, 7], [3, 8], [7, 3], [4, 3], [8, 4], [3, 8], [3, 4], [8, 4],
    #              [3, 8], [4, 3], [3, 9], [4, 3], [9, 4], [3, 9], [3, 4], [9, 4], [3, 9], [4, 3], [4, 9], [9, 0], [8, 9],
    #              [0, 8], [0, 9], [8, 9], [0, 8], [9, 0], [8, 0], [9, 8], [0, 9], [0, 8], [9, 8], [0, 9], [8, 0]]
    # input_filename = 'E:/python/Distributed_quantum_circuits_scheduling/qasm/8bitadder2.qasm'
    # gate_list = list_str_to_int(converter_circ_from_qasm(input_filename))
    # gate_list =[[15, 14], [15, 13], [15, 12], [15, 11], [15, 10], [15, 9], [15, 8], [15, 7], [15, 6], [15, 5], [15, 4], [15, 3], [15, 2], [15, 1], [15, 0], [14, 13], [14, 12], [14, 11], [14, 10], [14, 9], [14, 8], [14, 7], [14, 6], [14, 5], [14, 4], [14, 3], [14, 2], [14, 1], [14, 0], [13, 12], [13, 11], [13, 10], [13, 9], [13, 8], [13, 7], [13, 6], [13, 5], [13, 4], [13, 3], [13, 2], [13, 1], [13, 0], [12, 11], [12, 10], [12, 9], [12, 8], [12, 7], [12, 6], [12, 5], [12, 4], [12, 3], [12, 2], [12, 1], [12, 0], [11, 10], [11, 9], [11, 8], [11, 7], [11, 6], [11, 5], [11, 4], [11, 3], [11, 2], [11, 1], [11, 0], [10, 9], [10, 8], [10, 7], [10, 6], [10, 5], [10, 4], [10, 3], [10, 2], [10, 1], [10, 0], [9, 8], [9, 7], [9, 6], [9, 5], [9, 4], [9, 3], [9, 2], [9, 1], [9, 0], [8, 7], [8, 6], [8, 5], [8, 4], [8, 3], [8, 2], [8, 1], [8, 0], [7, 6], [7, 5], [7, 4], [7, 3], [7, 2], [7, 1], [7, 0], [6, 5], [6, 4], [6, 3], [6, 2], [6, 1], [6, 0], [5, 4], [5, 3], [5, 2], [5, 1], [5, 0], [4, 3], [4, 2], [4, 1], [4, 0], [3, 2], [3, 1], [3, 0], [2, 1], [2, 0], [1, 0]]
    gate_list = [[2, 1], [6, 2], [1, 6], [1, 2], [6, 2], [1, 6], [2, 1], [6, 5], [3, 6], [5, 3], [5, 6], [3, 6], [5, 3], [6, 5],
     [2, 1], [4, 2], [1, 4], [1, 2], [4, 2], [1, 4], [2, 1], [2, 4], [4, 0], [3, 4], [0, 3], [0, 4], [3, 4], [0, 3],
     [4, 0], [3, 0], [4, 3], [0, 4], [0, 3], [4, 3], [0, 4], [3, 0], [1, 7], [2, 7], [5, 9], [6, 9], [7, 5], [9, 7],
     [5, 9], [5, 7], [9, 7], [5, 9], [7, 5], [6, 5], [9, 6], [5, 9], [5, 6], [9, 6], [5, 9], [6, 5], [5, 8], [6, 5],
     [8, 6], [5, 8], [5, 6], [8, 6], [5, 8], [6, 5], [6, 8], [8, 0], [9, 8], [0, 9], [0, 8], [9, 8], [0, 9], [8, 0],
     [9, 0], [8, 9], [0, 8], [0, 9], [8, 9], [0, 8], [9, 0]]

    print(len(gate_list))

    circuit_qubit = max([max(row) for row in gate_list]) +1

    new_circuit(gate_list,circuit_qubit)

    # new_qft_circuit(4)

    img = Image.open(img_path)
    img.show()