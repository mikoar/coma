import pandas as pd
import argparse
from read_files import read_alignments_file, read_all_files, read_segments_file
import re
from write_indel_files import write_indel_file
from src.args import Args
from src.extensions.extension import Extension
from src.extensions.messages import AlignmentResultRowMessage, MultipleAlignmentResultRowsMessage
from src.program import Program


def main():
    parser = argparse.ArgumentParser(description="Looks for indels in alignment file based on segments conflicts")
    parser.add_argument("-r", "--reference", dest="referenceFile", type=str, required=True)
    parser.add_argument("-q", "--query", dest="queryFile", type=str, required=True)
    parser.add_argument("-a", "--alignedOutput", dest="alignedFile", type=str, required=True,
                        help="Name of output file with COMA alignments")
    parser.add_argument("-s", "--segmentsOutput", dest="segmentsFile", type=str, required=True,
                        help="Name of output file with multiple segments created during COMA workflow")
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?",
                        type=str, default="segment_indels.txt")
    parser.add_argument("-ma", "--secondaryMargin", dest="secondaryMargin", type=int, default=16000,
                        help="The number of base pairs by which the peak from initial cross-correlation "
                        "seeding is extended in both directions to serve as an input for the "
                        "second cross-correlation run.")
    parser.add_argument("-pt", "--peakHeightThreshold", dest="peakHeightThreshold", type=float, default=27,
                        help="Minimum second cross-correlation peak height to qualify for aligned pairs search.")
    parser.add_argument("-er", "--endReachingScore", dest="endReachingScore", type=int, default=0,
                            help="Score value given to an alignment reaching query start or end.")
    parser.add_argument("-sc", "--segmentCombinePenalty", dest="segmentCombinePenalty", type=int, default=0,
                        help="Constant component of penalty for combining segment.")
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
                        help="Minimum score of a segment/alignment.")
    parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=1200,
                        help="Alignment segments can be split into two if their score drops below this threshold.")
    parser.add_argument("-diff", "--maxDifference", dest="maxDifference", type=int, default=100000,
                        help="Multiple alignments of the same query will be joined if difference between their "
                                "reference positions is less or equal this parameter.")
    parser.add_argument("-pb", "--disableProgressBar", dest="disableProgressBar", action="store_true",
                        help="Disables the progress bar.")

    args = parser.parse_args()  # type: ignore
    run(args)



class SegmentsCatcher(Extension):
    messageType = AlignmentResultRowMessage

    def __init__(self, filePath: str):
        self.filePath = filePath

    def handle(self, message: AlignmentResultRowMessage):
        peak = message.correlation.maxPeak
        if peak:
            with open(self.filePath, "a") as f:
                f.write(
                        f"{message.correlation.query.moleculeId};"
                        f"{message.correlation.reverseStrand};"
                        f"{len(message.alignment.segments)};"
                        f"{message.alignment.segments};"
                        f"{message.alignment.queryStartPosition};"
                        f"{message.alignment.queryEndPosition};"
                        f"{message.alignment.referenceStartPosition};"
                        f"{message.alignment.referenceEndPosition};"
                        f"{message.alignment.confidence};"
                        f"{peak.score:.2f};"
                        f"{peak.position}\n")
                
class MultipleSegmentsCatcher(Extension):
    messageType = MultipleAlignmentResultRowsMessage

    def __init__(self, filePath: str):
        self.filePath = filePath

    def handle(self, message: AlignmentResultRowMessage):
        peak = message.correlation.maxPeak
        if peak:
            with open(self.filePath, "a") as f:
                f.write(
                        f"{message.correlation.query.moleculeId};"
                        f"{message.correlation.reverseStrand};"
                        f"{len(message.alignment.segments)};"
                        f"{message.alignment.segments};"
                        f"{message.alignment.queryStartPosition};"
                        f"{message.alignment.queryEndPosition};"
                        f"{message.alignment.referenceStartPosition};"
                        f"{message.alignment.referenceEndPosition};"
                        f"{peak.score:.2f};"
                        f"{peak.position}\n")

def find_conflict_place(multiple_segments: dict, alignmentFile: str) -> dict:
    """Function used to find places where segments where joined

    :param multiple_segments: Dictionary with molecules which had multiple segments
    :type multiple_segments: dict
    :param alignmentFile: COMA aligned file
    :type alignmentFile: str
    :return: Dictionary with identified places where segments where joined
    :rtype: dict
    """
    breakage_places = {}
    alignments = read_alignments_file(alignmentFile)
    for chromosme in alignments.values():
        for alignment in chromosme:
            # Analyze only alignments with multiple segments
            if alignment.queryId in multiple_segments.keys():
                added = False
                if len(multiple_segments[alignment.queryId]) == 2:
                    if str(alignment.alignedPairs[0]) == multiple_segments[alignment.queryId][0][0]:
                        first = multiple_segments[alignment.queryId][0]
                    elif str(alignment.alignedPairs[0]) == multiple_segments[alignment.queryId][1][0]:
                        first = multiple_segments[alignment.queryId][1]
                    else:
                        added = True
                    if added == False:
                        for index, pair in enumerate(first):
                            if str(alignment.alignedPairs[index]) != pair and added is False:
                                breakage_places[int(alignment.queryId)] = [[index, pair]]
                                added = True
                    if added == False:
                        breakage_places[int(alignment.queryId)] = [[index, pair]]
                elif len(multiple_segments[alignment.queryId]) > 2:
                    # Look for matching alignments
                    for segments_index in range((len( multiple_segments[alignment.queryId])-1)):
                        if str(alignment.alignedPairs[0]) == multiple_segments[alignment.queryId][segments_index][0]:
                            first = multiple_segments[alignment.queryId][segments_index]
                            if added == False:
                                for index, pair in enumerate(first):
                                    if str(alignment.alignedPairs[index]) != pair and added is False:
                                        breakage_places[int(alignment.queryId)] = [[index, pair]]
                                        added = True
                            # No breakage found during iteration of aligned pairs and found segment
                            if added == False:
                                breakage_places[int(alignment.queryId)] = [[index, pair]]
                                # Looking for other breakage points- make sure that alignment extends over first one
                                if len(alignment.alignedPairs) < index +1:
                                    if str(alignment.alignedPairs[index + 1]) == multiple_segments[alignment.queryId][segments_index+1][0]:
                                        added = False
                                        for second_index, pair in enumerate(multiple_segments[alignment.queryId][segments_index+1]):
                                            if str(alignment.alignedPairs[index + second_index + 1]) != pair and added is False:
                                                breakage_places[int(alignment.queryId)].append([index + second_index + 1, pair])
                                                added = True
                                        if added == False:
                                            breakage_places[int(alignment.queryId)].append([index + second_index + 1, pair])
    return breakage_places

def find_segments_indels(small_df: pd.DataFrame) -> dict:
    """Function used to create dictionary with multiple
    segments for molecules which had them

    :param small_df: DataFrame with observed segments for molecules
    :type small_df: pd.DataFrame
    :return: Dictionary with multiple segments for molecules which had them
    :rtype: dict
    """
    segment_d = {}
    for _, row in small_df.iterrows():
        segments = row["segment"]
        segments = segments.split("],")
        segments = [i for i in segments if i not in [' score: 0.0, positions: [', ' score: 0.0, positions: []]']]
        if len(segments) > 1:
            for i in segments:
                if row["queryId"] not in segment_d.keys():
                    if len([pair for pair in re.findall(r'\(.*?\)', i) if "-" not in pair]) > 0:
                        segment_d[row["queryId"]] = [[pair for pair in re.findall(r'\(.*?\)', i) if "-" not in pair]]
                else:
                    if len([pair for pair in re.findall(r'\(.*?\)', i) if "-" not in pair]) > 0:
                        segment_d[row["queryId"]].extend([[pair for pair in re.findall(r'\(.*?\)', i) if "-" not in pair]])
    return segment_d

def look_for_indels_in_breakage(alignment_dict: dict, reference_dict: dict,
                                query_dict: dict, breakage_dict: dict) -> dict:
    """Function ued to identify indels in places where segments are joined

    :param alignment_dict: Dictionary with aligned molecules
    :type alignment_dict: dict
    :param reference_dict: Dictionary with reference molecules
    :type reference_dict: dict
    :param query_dict: Dictionary with query molecules
    :type query_dict: dict
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
            if q_id in breakage_dict.keys():
                for i in breakage_dict[q_id]:
                    breakage_place = i
                    if len(alignment.alignedPairs) > breakage_place[0] + 1:
                        next_pair = alignment.alignedPairs[breakage_place[0] + 1]
                        breakage_pair = alignment.alignedPairs[breakage_place[0]]

                        r_label_s = reference_dict[r_id].positions[breakage_pair.reference.siteId - 1]
                        q_label_s = query_dict[q_id].positions[breakage_pair.query.siteId - 1]
                        r_label_e = reference_dict[r_id].positions[next_pair.reference.siteId - 1]
                        q_label_e = query_dict[q_id].positions[next_pair.query.siteId - 1]

                        diff = abs(r_label_s - r_label_e) - abs(q_label_s - q_label_e)
                        if abs(diff) > 100 and abs(diff) < 100000:
                            if diff < -100:
                                indels["insertion"].append(
                                        ["insertion", alignment.referenceId, r_label_s, r_label_e,
                                        q_id, q_label_s, q_label_e, diff])
                            else:
                                indels["deletion"].append(
                                        ["deletion", alignment.referenceId, r_label_s, r_label_e,
                                        q_id, q_label_s, q_label_e, diff])
    return indels



def run(args):
    segments_file = args.segmentsFile
    with open(segments_file, "w") as f:
        f.write("queryId;reverseStrand;segmentNb;segment;segmentQueryStart;segmentQueryEnd;segmentReferenceStart;segmentReferenceEnd;alignmentConfidence;score;peakPosition\n")

    action_args = []
    if args.diagnosticsEnabled:
        action_args.append("-D")
    if args.disableProgressBar:
        action_args.append("-pb")
    program_args = Args.parse([
        "-q", args.queryFile,
        "-r", args.referenceFile,
        "-o", args.alignedFile,
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
        "-er", str(args.endReachingScore),
        "-sc", str(args.segmentCombinePenalty),
        "-sj", str(args.segmentJoinMultiplier),
        "-ss", str(args.sequentialityScore),
    ] + action_args)

    coma = Program(program_args, [SegmentsCatcher(segments_file)])
    coma.run()

    small_df = read_segments_file(segments_file)
    more_segments = find_segments_indels(small_df)
    
    with open(segments_file, "w") as f:
        f.write("queryId;reverseStrand;segmentNb;segment;segmentQueryStart;segmentQueryEnd;segmentReferenceStart;segmentReferenceEnd;alignmentConfidence;score;peakPosition\n")
    small_df.to_csv(segments_file, index=False)

    breakage_points = find_conflict_place(more_segments, args.alignedFile)

    reference_dict, query_dict, alignment_referenceId = read_all_files(reference_file=args.referenceFile,
                                                                       alignment_file=args.alignedFile,
                                                                       query_file=args.queryFile)

    indels_dict = look_for_indels_in_breakage(alignment_referenceId, reference_dict,
                                              query_dict, breakage_points)

    write_indel_file(indels_dict, args.alignedFile, file_name=args.outputFile)

if __name__ == '__main__':
    main()
