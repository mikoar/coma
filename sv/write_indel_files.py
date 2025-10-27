import operator

def cluster_indels(list_indels:list, blur:int=30000) -> list:
    """Identify and count same indels found in different query molecules

    :param list_indels: list onf indels of one type
    :type list_indels: list
    :param blur: PositionBlur which should be used during clustering, defaults to 30000
    :type blur: int, optional
    :return: Clustered lines
    :rtype: list
    """
    new_list = []
    if len(list_indels) > 0:
        first_line = list_indels[0]
        new_list = [first_line + [1]]
        for line in list_indels[1:]:
            if abs(line[3] - new_list[-1][3]) <= blur:
                if line[0:2] == new_list[-1][0:2]:
                    prev_line = new_list[-1]
                    if len(prev_line) == 8:
                        prev_line.append(1)
                    else:
                        prev_line[8] = prev_line[8] + 1
                        prev_line[2] = min(prev_line[2], line[2])
                        prev_line[3] = max(prev_line[3], line[3])
                        prev_line[7] = (prev_line[7]+ line[7])/2
                    prev_line[4] = str(prev_line[4]) + "," + str(line[4])
                    new_list[-1] = prev_line
            elif abs(line[2] - new_list[-1][2]) <= blur and abs(line[3] - new_list[-1][3]) <= blur:
                if line[0:2] == new_list[-1][0:2]:
                    prev_line = new_list[-1]
                    if len(prev_line) == 8:
                        prev_line.append(1)
                    else:
                        prev_line[8] = prev_line[8] + 1
                        prev_line[2] = min(prev_line[2], line[2])
                        prev_line[3] = max(prev_line[3], line[3])
                        prev_line[7] = (prev_line[7]+ line[7])/2
                    prev_line[4] = str(prev_line[4]) + "," + str(line[4])
                    new_list[-1] = prev_line
            else:
                new_list.append(line + [1])
    return new_list

def write_indel_file(indels_dict:dict, alignment_file_name:str, file_name:str="indels.txt"):
    """Function used to write indels files

    :param indels_dict: Dictionary with identified indels
    :type indels_dict: dict
    :param alignment_file_name: Name of alignments XMAP file
    :type alignment_file_name: str
    :param file_name: Name of output file, defaults to "indels.txt"
    :type file_name: str, optional
    """
    lines_deletions = list(indels_dict.values())[1]
    lines_insertions = list(indels_dict.values())[0]
    lines_sorted_del = sorted(lines_deletions, key=operator.itemgetter(1, 3))
    lines_sorted_ins = sorted(lines_insertions, key=operator.itemgetter(1, 3))

    lines_deletions_sorted = cluster_indels(lines_sorted_del, blur=-1)
    lines_insertions_sorted = cluster_indels(lines_sorted_ins, blur=-1)
    lines_sorted = sorted(lines_deletions_sorted + lines_insertions_sorted,
                          key=operator.itemgetter(1, 3))

    with open(file_name, "w") as f:
        f.write("#" + alignment_file_name + "\n")
        f.write("#Type \t Chromosome \t RefStart \t RefStop \t QueryId \t QueryStart \t QueryStop \t Length \t Count \n")
        for line in lines_sorted:
            line = [str(i) for i in line]
            line = "\t".join(line)
            f.write(line + "\n")
