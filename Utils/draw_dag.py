from qiskit import QuantumRegister, ClassicalRegister, AncillaRegister, QuantumCircuit, dagcircuit
from qiskit.dagcircuit import DAGCircuit
from qiskit.converters import circuit_to_dag
from qiskit.tools.visualization import dag_drawer

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
            circuit.x(gate_list[i][0])
    return circuit

if __name__ == '__main__':
    # gate_list = [[1, 2, 4], [3, 4, 5], [1, 2, 6], [2, 6], [0, 6, 5], [0, 5, 6], [1, 7], [2, 7], [3, 8], [4, 8],
    #              [3, 7, 8], [3, 4, 8], [3, 9], [3, 4, 9], [4, 9],  [0, 9, 8], [0, 8, 9]]
    gate_list = [[0,3],[2,3],[3,2],[0,3],[0,2],[3,2],[2,3],[1,2],[4,2],[2,4],[1,2],[1,4],[2,4],[4,2],[1,3],[5,3],[3,5],[1,3],[1,5],[3,5],[5,3]]
    circuit_qubit = max([max(row) for row in gate_list]) + 1
    circuit = new_circuit(gate_list, circuit_qubit)
    print(circuit)
    # circuit = QuantumCircuit(5)
    # circuit.ccx(0, 4, 2)
    # circuit.cx(2, 3)
    # circuit.cx(2, 4)
    # circuit.cx(1, 2)

    dag = circuit_to_dag(circuit)
    print(type(dag))
    # dag_drawer(dag)
