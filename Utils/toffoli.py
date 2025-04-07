
'''toffoli门分解
   分解成V门
'''
def toffoli_decompose_v(gate_list):
    decompose_gate_list = []
    for i in range(len(gate_list)):
        if len(gate_list[i]) < 3: # 单量子门或者cnot门
            decompose_gate_list.append(gate_list[i])
        # toffoli 的两种分解方法
        if len(gate_list[i]) == 3: # toffoli门
            decompose_gate_list.append([gate_list[i][0], gate_list[i][2]])
            decompose_gate_list.append([gate_list[i][0], gate_list[i][1]])
            decompose_gate_list.append([gate_list[i][1], gate_list[i][2]])
            decompose_gate_list.append([gate_list[i][0], gate_list[i][1]])
            decompose_gate_list.append([gate_list[i][1], gate_list[i][2]])
        # 两种分解方式
        # if len(gate_list[i]) == 3: # toffoli门
        #     decompose_gate_list.append([gate_list[i][1], gate_list[i][2]])
        #     decompose_gate_list.append([gate_list[i][0], gate_list[i][1]])
        #     decompose_gate_list.append([gate_list[i][1], gate_list[i][2]])
        #     decompose_gate_list.append([gate_list[i][0], gate_list[i][1]])
        #     decompose_gate_list.append([gate_list[i][0], gate_list[i][2]])
    # print(decompose_gate_list)
    return decompose_gate_list



'''toffoli门分解
   分解成cliford门
   只包含双量子门
'''
def toffoli_decompose_c(gate_list):
    decompose_gate_list = []
    for i in range(len(gate_list)):
        if len(gate_list[i]) < 3: # 单量子门或者cnot门
            decompose_gate_list.append(gate_list[i])
        if len(gate_list[i]) == 3: # toffoli门
            decompose_gate_list.append([gate_list[i][1], gate_list[i][2]])
            decompose_gate_list.append([gate_list[i][0], gate_list[i][2]])
            decompose_gate_list.append([gate_list[i][1], gate_list[i][2]])
            decompose_gate_list.append([gate_list[i][0], gate_list[i][2]])
            decompose_gate_list.append([gate_list[i][0], gate_list[i][1]])
            decompose_gate_list.append([gate_list[i][0], gate_list[i][1]])
    return decompose_gate_list