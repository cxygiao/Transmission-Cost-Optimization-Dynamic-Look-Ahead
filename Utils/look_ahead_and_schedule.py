from Utils.division_by_irregularity import split_into_k_parts, transfer_qubit_list_by_gate_list_based_look_ahead, \
    transfer_qubit_list_by_gate_list_based_look_ahead2
from Utils.generate_partitions import random_line
from Utils.min_global_gate_num import Initial_line_sequence, k_change_gate_list_by_line_sequence
from 新建文件夹.min_st_num import statistics_gate_labels
from Utils.read_qasm import get_data
from Utils.readreals import real_to_cnot_list
from Utils.transmission_cost_calculation_more import count_transfer_queue, new_global_gate_list, \
    judge_is_same_partation, judge_middle_gates_is_affect_transfer_qubit, count_transfer_queue_by_consider_last

'''K分区下计算最小传输代价 前瞻'''
def k_count_min_st_num_ahead(gate_list,cut_list,line_sequence_list):
    initial_line_sequence = Initial_line_sequence(gate_list)
    ST = []
    ST_schedule = []
    lenl = len(line_sequence_list)
    print(lenl)
    mini = len(gate_list)*2
    mini_schedule = mini
    for i in range(lenl):
        print(i, end=':')
        print(line_sequence_list[i],end='')
        new_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[i])
        st,st_schedule = direct_calculation_of_tc_look_ahead_and_schenule(new_gate_list, cut_list)
        ST.append(st)
        ST_schedule.append(st_schedule)
        print('当前传输代价最小：',mini,' ','当前传输代价最小：',mini_schedule)
        if i%100 == 0:
            mini = min(ST)
            mini_schedule = min(ST_schedule)
    print('最小传输代价（未调度）：', min(ST))
    print('最小传输代价（调度）：', min(ST_schedule))
    # min_st_gate_list = k_change_gate_list_by_line_sequence(initial_line_sequence,gate_list,line_sequence_list[ST.index(min(ST))])
    # 绘制电路
    # draw_color_circuit(min_st_gate_list)
    return ST

'''初始解直接搜索 通过队列'''
def direct_calculation_of_tc_look_ahead_and_schenule(gate_list,cut_list):
    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)

    global_gate_list = new_global_gate_list(statistics_gate_labels_list)
    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list0 = transfer_qubit_list_by_gate_list_based_look_ahead(gate_list,
                                                                                     statistics_gate_labels_list)
    initial_transfer_qubit_list1 = transfer_qubit_list_by_gate_list_based_look_ahead2(gate_list,
                                                                                      statistics_gate_labels_list)
    # 统计合并传输队列
    # print('######################################################################################################')
    initial_transfer_queue0 = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list0,
                                                   statistics_gate_labels_list)
    initial_transfer_queue1 = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list1,
                                                   statistics_gate_labels_list)
    st0 = len(initial_transfer_queue0) * 2
    st1 = len(initial_transfer_queue1) * 2

    print('初始传输代价：' + str(min(st0, st1)), end=' ')

    # print('######################################################################################################')
    st_schedule0 = schedule(initial_transfer_queue0,initial_transfer_qubit_list0,global_gate_list)
    st_schedule1 = schedule(initial_transfer_queue1, initial_transfer_qubit_list1, global_gate_list)

    print('传输调度后传输代价：' + str(min(st_schedule0, st_schedule1)),end=' ')
    return min(st0, st1),min(st_schedule0, st_schedule1)


#调度
def schedule(initial_transfer_queue,initial_transfer_qubit_list,global_gate_list):

    st = len(initial_transfer_queue) * 2
    # 如果两个队列间隔的一些量子门与传输量子位无相同量子位，且量子位传输列表位于同一个分区，则回传和正传可省略
    for i in range(len(initial_transfer_queue) - 1):

        # 量子位传输列表位于同一个分区
        gate_i = initial_transfer_queue[i][0]  # g2
        gate_num_i = global_gate_list.index(gate_i)  # ['g2', 'g6', 'g9'] -> 0
        transfer_qubit_i = int(get_data(initial_transfer_qubit_list[gate_num_i])[0])  #
        gate_j = initial_transfer_queue[i + 1][0]  # g6
        gate_num_j = global_gate_list.index(gate_j)
        transfer_qubit_j = int(get_data(initial_transfer_qubit_list[gate_num_j])[0])

        gate_i_next_name = initial_transfer_queue[i][-1]  # ['g2', 'g6', 'g9'] -> g9
        gate_i_next = int(gate_i_next_name[
                          1:]) + 1  # [[2, 0], [0, 2], [2, 5], [3, 4], [1, 2], [2, 1], [1, 4], [4, 3], [0, 2], [4, 1]] -> 3
        gate_j_front = int(gate_j[1:])  # -> 6
        middle_gate_list = gate_list[gate_i_next:gate_j_front]  # 左闭右开，所以 gate_j_front 不用 +1

        # 两个队列间隔的一些量子门与传输量子位无相同量子位   起始点：第一个队列中第一个门 结束点：第二个队列中第一个门
        if judge_is_same_partation(transfer_qubit_i, transfer_qubit_j,
                                   cut_list) == 1 and judge_middle_gates_is_affect_transfer_qubit(middle_gate_list,
                                                                                       transfer_qubit_i) == 0:
            # print(gate_i)
            # print(gate_j)
            st -= 1

    return st


'''初始解直接搜索 通过队列'''
def direct_calculation_of_tc_look_ahead_new(gate_list,cut_list):
    # 全局门标签
    statistics_gate_labels_list = statistics_gate_labels(gate_list, cut_list)

    global_gate_list = new_global_gate_list(statistics_gate_labels_list)
    # 量子位传输列表 ['q0', 'q0', 'q2', 'q0', 'q4', 'q1', 'q0', 'q4']
    initial_transfer_qubit_list0 = transfer_qubit_list_by_gate_list_based_look_ahead(gate_list,
                                                                                     statistics_gate_labels_list)
    initial_transfer_qubit_list1 = transfer_qubit_list_by_gate_list_based_look_ahead2(gate_list,
                                                                                      statistics_gate_labels_list)
    # 统计合并传输队列
    # print('######################################################################################################')
    # initial_transfer_queue0 = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list0,
    #                                                statistics_gate_labels_list)
    # initial_transfer_queue1 = count_transfer_queue(gate_list, cut_list, initial_transfer_qubit_list1,
    #                                                statistics_gate_labels_list)
    # st0 = len(initial_transfer_queue0) * 2
    # st1 = len(initial_transfer_queue1) * 2
    # 考虑最后的量子门的回传
    st0 = count_transfer_queue_by_consider_last(gate_list, cut_list, initial_transfer_qubit_list0,
                                                   statistics_gate_labels_list)
    st1 = count_transfer_queue_by_consider_last(gate_list, cut_list, initial_transfer_qubit_list1,
                                                   statistics_gate_labels_list)
    # print(initial_transfer_queue0)
    # print('初始传输代价：' + str(min(st0, st1)), end=' ')
    return min(st0, st1)

if __name__ == '__main__':
    input_filename = '../reals/new_real/8bitadder.real'
    gate_list = real_to_cnot_list(input_filename)
    # gate_list = qasm_to_cnot_list(input_filename)
    qubit = max([max(row) for row in gate_list]) + 1

    cut_list = split_into_k_parts(qubit, 6)
    # cut_list = [5,4,3]
    # 遍历排列组合
    # line_sequence_list = generate_line(cut_list)
    line_sequence_list = random_line(cut_list, 1000000)
    st = k_count_min_st_num_ahead(gate_list,cut_list,line_sequence_list)

    print(input_filename)
    print(cut_list)