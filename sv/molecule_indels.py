from read_files import read_alignments_file, read_all_files
from write_indel_files import write_indel_file
import argparse
from src.parsers.xmap_reader import BionanoAlignment


def main():
    parser = argparse.ArgumentParser(description="Looks for indels in joined alignment file")
    parser.add_argument("-r", "--reference", dest="referenceFile", type=str, required=True)
    parser.add_argument("-q", "--query", dest="queryFile", type=str, required=True)
    parser.add_argument("-j", "--joined", dest="joinedFile", type=str, required=True,
                        help="File with joined alignments")
    parser.add_argument("-f", "--first", dest="firstFile", type=str, required=True,
                        help="File with FIRST pass alignments")
    parser.add_argument("-s", "--second", dest="secondFile", type=str, required=True,
                        help="File with SECOND pass alignments")
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=str, default="joined_indels.txt")

    args = parser.parse_args()  # type: ignore
    run(args)


def get_joined_ids(joined_file_name: str) -> set:
    """Function used to get ides of the joined molecules

    :param joined_file_name: Name of the joined aligned XMAP file
    :type joined_file_name: str
    :return: Set with molecules ids present in XMAP file
    :rtype: set
    """
    joined_query_ids = set()
    with open(joined_file_name, 'r') as f:
        lines = f.readlines()
    lines = [line.rstrip() for line in lines if not line.startswith("#")]
    for line in lines:
        line = line.split("\t")
        joined_query_ids.add(int(line[1]))
    return joined_query_ids

def add_alignment_to_dict(current_dict: dict, alignment: BionanoAlignment, alignment_type: str) -> dict:
    """Function used to update alignment dictionary

    :param current_dict: Alignment dictionary
    :type current_dict: dict
    :param alignment: Alignment which will be added to the dictionary
    :type alignment: BionanoAlignment
    :param alignment_type: Whether the alignments was obtained during FIST or SECOND PASS
    :type alignment_type: str
    :return: Updated dictionary
    :rtype: dict
    """
    current_dict[alignment.queryId][alignment_type] = alignment
    return current_dict


def get_separate_alignments(original_alignment_file: str, rests_alignment_file: str,
                            selected_ids: set) -> dict:
    """Function used to create dictionary with
    FIRST and SECOND PASS alignments for each molecule

    :param original_alignment_file: Name of FIRST PASS file
    :type original_alignment_file: str
    :param rests_alignment_file: Name of SECOND PASS file
    :type rests_alignment_file: str
    :param selected_ids: List of ids of joined molecules
    :type selected_ids: set
    :return: Dictionary with FIRS and SECOND PASS alignments
    for each molecule present in the joined file
    :rtype: dict
    """
    selected_ids_alignments = {i : {} for i in selected_ids}
    original_alignments = read_alignments_file(original_alignment_file)
    rests_alignments = read_alignments_file(rests_alignment_file)
    for chromosome, values in original_alignments.items():
        for value in values:
            if value.queryId in selected_ids:
                selected_ids_alignments = add_alignment_to_dict(selected_ids_alignments,
                                                                value, "original")
        for value in rests_alignments[chromosome]:
            if value.queryId in selected_ids:
                selected_ids_alignments = add_alignment_to_dict(selected_ids_alignments,
                                                                value, "rest")
    return selected_ids_alignments


def find_conflict_place(joint_alignments: str, multiple_alignments: dict) -> dict:
    """Function used to identify breakage point in joined alignments

    :param joint_alignments: Name of the XMAP file with joined alignments
    :type joint_alignments: str
    :param multiple_alignments: Dictionary with FIRST and SECOND pass alignments
    :type multiple_alignments: dict
    :return: Dictionary with identifies breakage places for all molecules
    :rtype: dict
    """

    breakage_places = {}
    joined = read_alignments_file(joint_alignments)
    for chromosme in joined.values():
        for alignment in chromosme:
            added = False
            if alignment.alignedPairs[0] == multiple_alignments[alignment.queryId]["original"].alignedPairs[0]:
                first = multiple_alignments[alignment.queryId]["original"]
            else:
                first = multiple_alignments[alignment.queryId]["rest"]
            for index, pair in enumerate(first.alignedPairs):
                if alignment.alignedPairs[index] != pair and added is False:
                    breakage_places[int(alignment.queryId)] = [index, pair]
                    added = True
                    break
            if added == False:
                breakage_places[int(alignment.queryId)] = [index, pair]
    return breakage_places


def look_for_indels_in_breakage(alignment_dict: dict, r_dict: dict,
                                q_dict: dict, breakage_dict: dict) -> dict:
    """Function ued to identify indels in places where molecules are joined

    :param alignment_dict: Dictionary with joined alignments
    :type alignment_dict: dict
    :param r_dict: Dictionary with reference molecules
    :type r_dict: dict
    :param q_dict: Dictionary with query molecules
    :type q_dict: dict
    :param breakage_dict: Dictionary with information about the breakage points
    :type breakage_dict: dict
    :return: Dictionary with identified insertions and deletions
    :rtype: dict
    """
    indels = {"insertion" : [],
              "deletion" : []}
    for chromosome_alignments in alignment_dict.values():
        for alignment in chromosome_alignments:
            q_id = alignment.queryId
            r_id = alignment.referenceId

            breakage_place = breakage_dict[q_id]
            next_pair = alignment.alignedPairs[breakage_place[0] + 1]

            r_label_s = r_dict[r_id].positions[breakage_place[1].reference.siteId - 1]
            q_label_s = q_dict[q_id].positions[breakage_place[1].query.siteId - 1]
            r_label_e = r_dict[r_id].positions[next_pair.reference.siteId - 1]
            q_label_e = q_dict[q_id].positions[next_pair.query.siteId - 1]

            diff = abs(r_label_s - r_label_e) - abs(q_label_s - q_label_e)
            if abs(diff) > 2000 and abs(diff) < 100000:
                if diff < -2000:
                    indels["insertion"].append(
                            ["insertion", alignment.referenceId, r_label_s, r_label_e,
                            q_id, q_label_s, q_label_e, diff])
                else:
                    indels["deletion"].append(
                            ["deletion", alignment.referenceId, r_label_s, r_label_e,
                            q_id, q_label_s, q_label_e, diff])
    return indels


def run(args):
    joined_queries_id = get_joined_ids(args.joinedFile)
    queryID_alignments_dict = get_separate_alignments(args.firstFile, args.secondFile,
                                                      joined_queries_id)
    breakage_points = find_conflict_place(args.joinedFile, queryID_alignments_dict)

    reference_dict, query_dict, alignment_referenceId = read_all_files(reference_file=args.referenceFile,
                                                                       alignment_file=args.joinedFile,
                                                                       query_file=args.queryFile)
    indels_dict = look_for_indels_in_breakage(alignment_referenceId, reference_dict,
                                              query_dict, breakage_points)

    write_indel_file(indels_dict, args.joinedFile, file_name=args.outputFile)

if __name__ == '__main__':
    main()
