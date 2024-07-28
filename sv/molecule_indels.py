from read_files import read_alignments_file, read_all_files
from write_indel_files import write_indel_file
import argparse
from src.parsers.xmap_reader import BionanoAlignment
from src.args import Args
from src.program import Program


def main():
    parser = argparse.ArgumentParser(description="Looks for indels in joined alignment file created during included coma run")
    parser.add_argument("-r", "--reference", dest="referenceFile", type=str, required=True)
    parser.add_argument("-q", "--query", dest="queryFile", type=str, required=True)
    parser.add_argument("-j", "--joined", dest="joinedFile", type=str, default='aligned_all.xmap',
                        help="File with joined alignments")
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=str, default="joined_indels.txt")
    parser.add_argument("-ma", "--secondaryMargin", dest="secondaryMargin", type=int, default=16000,
                        help="The number of base pairs by which the peak from initial cross-correlation "
                        "seeding is extended in both directions to serve as an input for the "
                        "second cross-correlation run.")
    parser.add_argument("-pt", "--peakHeightThreshold", dest="peakHeightThreshold", type=float, default=27,
                        help="Minimum second cross-correlation peak height to qualify for aligned pairs search.")
    parser.add_argument("-sj", "--segmentJoinMultiplier", dest="segmentJoinMultiplier", type=float, default=1,
                        help="Multiplier applied to segment sequentiality scores.")
    parser.add_argument("-ss", "--sequentialityScore", dest="sequentialityScore", type=int, default=0,
                        help="Segment sequentiality scoring function version.")
    parser.add_argument("-D", "--diagnostics", dest="diagnosticsEnabled", action="store_true",
                        help="Draws cross-correlation and alignment plots. When used, 'alignmentFile' parameter "
                             "is required")
    parser.add_argument("-r1", "--primaryResolution", dest="primaryResolution", type=int, default=1400,
                        help="Scaling factor used to reduce the size of the vectorized form of the optical map "
                                "in the initial cross-correlation seeding step.")
    parser.add_argument("-b1", "--primaryBlur", dest="primaryBlur", type=int, default=1,
                        help="Extends each label in the vectorized form of the optical map in both directions "
                                "by given number of positions in the initial cross-correlation seeding step "
                                "in order to increase the chance of overlap. "
                                "Final width of each label is equal to 2b + 1.")
    parser.add_argument("-p", "--peaksCount", dest="peaksCount", type=int, default=3,
                        help="Number of peaks found for each query molecule against all reference molecules in the "
                                "first cross-correlation run that are selected for further steps - the second "
                                "cross-correlation run and alignment creation. Then the alignment with the highest "
                                "score is returned, one alignment record per query molecule at most.")
    parser.add_argument("-md", "--minPeakDistance", dest="minPeakDistance", type=int, default=20000,
                        help="Minimum distance between peaks identified in the initial cross-correlation. "
                                "For more details see parameter distance of scipy.signal._peak_finding.find_peaks.")
    parser.add_argument("-r2", "--secondaryResolution", dest="secondaryResolution", type=int, default=100,
                        help="Scaling factor used to reduce the size of the vectorized form of the optical map "
                                "in the second cross-correlation run.")
    parser.add_argument("-b2", "--secondaryBlur", dest="secondaryBlur", type=int, default=4,
                        help="Extends each label in the vectorized form of the optical map in both directions "
                                "by given number of positions in the second cross-correlation run "
                                "in order to increase the chance of overlap. "
                                "Final width of each label is equal to 2b + 1.")
    parser.add_argument("-d", "--maxPairDistance", dest="maxPairDistance", type=int, default=1500,
                        help="Maximum distance between aligned pairs relatively to the cross-correlation lag.")
    parser.add_argument("-sp", "--perfectMatchScore", dest="perfectMatchScore", type=int, default=1000,
                        help="Score value given to an aligned pair with 0 distance between reference and query "
                                "positions.")
    parser.add_argument("-dp", "--distancePenaltyMultiplier", dest="distancePenaltyMultiplier", type=float,
                        default=1., help="Multiplier applied to the distance between reference and query positions "
                                            "of an aligned pair that reduces the pair's score.")
    parser.add_argument("-su", "--unmatchedPenalty", dest="unmatchedPenalty", type=int, default=-250,
                        help="Penalty to a segment score for each unpaired reference or query position.")
    parser.add_argument("-ms", "--minScore", dest="minScore", type=int, default=1000,
                        help="Minimum score of a segment.")
    parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=1200,
                        help="Alignment segments can be split into two if their score drops below this threshold.")
    parser.add_argument("-diff", "--maxDifference", dest="maxDifference", type=int, default=100000,
                        help="Multiple alignments of the same query will be joined if difference between their "
                                "reference positions is less or equal this parameter.")
    parser.add_argument("-pb", "--disableProgressBar", dest="disableProgressBar", action="store_true",
                        help="Disables the progress bar.")

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
            # Make sure that alignment extends over breakage
            if len(alignment.alignedPairs) > breakage_place[0] + 1:
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
    action_args = []
    if args.diagnosticsEnabled:
        action_args.append("-D")
    if args.disableProgressBar:
        action_args.append("-pb")
    program_args = Args.parse([
        "-q", args.queryFile,
        "-r", args.referenceFile,
        "-o", args.joinedFile,
        "-r1", str(args.primaryResolution),
        "-b1", str(args.primaryBlur),
        "-p", str(args.peaksCount),
        "-md", str(args.minPeakDistance),
        "-r2", str(args.secondaryResolution),
        "-b2", str(args.secondaryBlur),
        "-d", str(args.maxPairDistance),
        "-sp", str(args.perfectMatchScore),
        "-dp", str(args.distancePenaltyMultiplier),
        "-su", str(args.unmatchedPenalty),
        "-ms", str(args.minScore),
        "-bs", str(args.breakSegmentThreshold),
        "-diff", str(args.maxDifference),
        "-ma", str(args.secondaryMargin),
        "-pt", str(args.peakHeightThreshold),
        "-sj", str(args.segmentJoinMultiplier),
        "-ss", str(args.sequentialityScore),
        "-oM", "all"
    ] + action_args)
    coma = Program(program_args)
    coma.run()
    
    joined_queries_id = get_joined_ids(args.joinedFile)
    queryID_alignments_dict = get_separate_alignments(str(args.joinedFile).replace(".xmap", "_1.xmap"), str(args.joinedFile).replace(".xmap", "_2.xmap"),
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
