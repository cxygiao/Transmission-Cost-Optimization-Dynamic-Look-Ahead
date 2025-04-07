from itertools import combinations
import random

def generate_partitions(remaining_letters, partition_sizes, current_partition, all_partitions):
    if len(partition_sizes) == 0:
        if len(remaining_letters) == 0:
            all_partitions.append(current_partition.copy())
            # print(current_partition.copy())
        return

    current_size = partition_sizes[0]
    for comb in combinations(remaining_letters, current_size):
        next_partition = current_partition.copy()
        next_partition.append(list(comb))
        next_remaining_letters = [letter for letter in remaining_letters if letter not in comb]

        generate_partitions(next_remaining_letters, partition_sizes[1:], next_partition, all_partitions)

def generate_line(partition_sizes):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9','$','#']

    letters = str2[0:sum(partition_sizes)]

    all_partitions = []

    generate_partitions(letters, partition_sizes, [], all_partitions)

    new_all_partitions = []
    for partition in all_partitions:
        f_partition = []
        for i in range(len(partition_sizes)):
            f_partition = f_partition+partition[i]
            str_partition = ''.join(f_partition)
        new_all_partitions.append(str_partition)

    return new_all_partitions



# def generate_partitions(remaining_letters, partition_sizes, current_partition, all_partitions):
#     if len(partition_sizes) == 0:
#         if len(remaining_letters) == 0:
#             line = "".join(["".join(row) for row in current_partition.copy()])
#             all_partitions.append(line)
#         return
#
#     current_size = partition_sizes[0]
#     for comb in itertools.combinations(remaining_letters, current_size):
#         next_partition = current_partition.copy()
#         next_partition.append(list(comb))
#         next_remaining_letters = [letter for letter in remaining_letters if letter not in comb]
#
#         generate_partitions(next_remaining_letters, partition_sizes[1:], next_partition, all_partitions)
#
#
# def generate_line2(partition_sizes):
#     str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
#             'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
#             'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
#
#     letters = str2[0:sum(partition_sizes)]
#
#     all_partitions = []
#
#     generate_partitions(letters, partition_sizes, [], all_partitions)
#
#     return all_partitions

def random_line(partition_sizes,random_num):
    str2 = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U',
            'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p',
            'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9','$','#']
    letters = str2[0:sum(partition_sizes)]
    all_partitions = []
    for i in range(random_num):
        random.shuffle(letters)
        all_partitions.append(''.join(letters))
    return all_partitions



if __name__ == '__main__':
    partition_sizes = [10,10]
    # all_partitions = generate_line(partition_sizes)
    # print(all_partitions)

    random_num = 10000
    random_partitions = random_line(partition_sizes, random_num)
    print(random_partitions)


