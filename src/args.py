from __future__ import annotations

import argparse
import sys
from typing import NamedTuple, TextIO, List, Literal


class Args(NamedTuple):
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO
    outputIndels: str
    outputRests: str
    outputSmap: str
    primaryResolution: int
    primaryBlur: int
    secondaryResolution: int
    secondaryBlur: int
    secondaryMargin: int
    referenceIds: List[int]
    queryIds: List[int]
    numberOfCpus: int | None
    minPeakDistance: int
    maxPairDistance: int
    peakHeightThreshold: float
    perfectMatchScore: int
    distancePenaltyMultiplier: int
    unmatchedPenalty: int
    minScore: int
    breakSegmentThreshold: int
    maxDifference: int
    diagnosticsEnabled: bool
    benchmarkAlignmentFile: TextIO
    peaksCount: int
    disableProgressBar: bool
    outputMode: Literal["best", "separate", "joined", "all", "single", "split-alignments"]
    segmentCombinePenalty: int
    segmentJoinMultiplier: float
    sequentialityScore: int
    endReachingScore: int
    minSubsequentScore: int
    scalingRange: float

    @staticmethod
    def parse(args: List[str] = None) -> Args:
        parser = argparse.ArgumentParser(description="Optical map aligner.")

        parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"), required=True,
                            help="Reference optical map in CMAP format file path.")

        parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"), required=True,
                            help="Query optical map in CMAP format file path.")

        parser.add_argument("-rId", "--referenceIDs", dest="referenceIds", type=int, nargs="*",
                            help="CMapId(s) of reference molecules to be used. Takes all if omitted.")

        parser.add_argument("-qId", "--queryIDs", dest="queryIds", type=int, nargs="*",
                            help="CMapId(s) of query molecules to be used. Takes all if omitted.")

        parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"),
                            default=sys.stdout,
                            help="XMAP output file path. Stdout is used if omitted.")

        parser.add_argument("-id", "--outputIndels", dest="outputIndels", nargs="?", type=str, default=None,
                            help="Segment indel output file path.")

        parser.add_argument("-rt", "--outputRests", dest="outputRests", nargs="?", type=str, default=None,
                            help="Alignment rest output file path.")

        parser.add_argument("-sm", "--smap", dest="outputSmap", nargs="?", type=str, default=None,
                            help="SMAP output file path.")

        parser.add_argument("-oM", "--outputMode", dest="outputMode", type=str,
                            default="single", choices=["best", "separate", "joined", "all", "single", "split-alignments"],
                            help="Mode which should be used while running aligner and creating output alignment file(s). "
                                "There are several possible options: "
                                "    'single' - single-pass running mode (default); "
                                "    'separate' - creates two separate files for alignments from first and second pass; "
                                "    'joined'- joins alignments when it is possible and saves rest to separate file; "
                                "    'best' - includes joined alignments when possible and best alignment based on "
                                "        confidence when joined option is not available; "
                                "    'all'- creates 3 files, one with joint alignments, and two with all obtained alignments; "
                                "    'split-alignments' - creates single file for all multi-pass alignments.")

        parser.add_argument("-sr", "--scalingRange", dest="scalingRange", type=float, default=0,
                            help="The extent to which a molecule can be stretched or compressed relative to its length "
                                "(default: no stretching).")

        parser.add_argument("-st", "--scalingStep", dest="scalingStep", type=int, default=None,
                            help="Stretching step of query molecule scaling (defaul: equals to parameter primaryResolution). ")

        parser.add_argument("-r1", "--primaryResolution", dest="primaryResolution", type=int, default=2000,
                            help="Size of sequence intervals represented by a single item in the vectorized form of the optical map "
                                 "in the initial cross-correlation seeding step.")

        parser.add_argument("-b1", "--primaryBlur", dest="primaryBlur", type=int, default=1,
                            help="Extends each label in the vectorized form of the optical map in both directions "
                                 "by given number of positions in the initial cross-correlation seeding step "
                                 "in order to increase the chance of overlap. "
                                 "Final width of each label is equal to (2*b1+1)*r1 bp.")

        parser.add_argument("-p1", "--primaryPeaksCount", dest="primaryPeaksCount", type=int, default=3,
                            help="Number of peaks found for each query molecule against all reference molecules in the "
                                 "first cross-correlation run that are selected for further steps - the second "
                                 "cross-correlation run and alignment creation. ")

        parser.add_argument("-md", "--minPeakDistance", dest="minPeakDistance", type=int, default=20000,
                            help="Minimum distance between peaks identified in the initial cross-correlation. "
                                 "For more details see parameter distance of scipy.signal._peak_finding.find_peaks.")

        parser.add_argument("-r2", "--secondaryResolution", dest="secondaryResolution", type=int, default=50,
                            help="Size of sequence intervals represented by a single item in the vectorized form of the optical map "
                                 "in the second cross-correlation run.")

        parser.add_argument("-b2", "--secondaryBlur", dest="secondaryBlur", type=int, default=4,
                            help="Extends each label in the vectorized form of the optical map in both directions "
                                 "by given number of positions in the second cross-correlation run "
                                 "in order to increase the chance of overlap. "
                                 "Final width of each label is equal to (2*b2+1)*r2 bp.")

        parser.add_argument("-p2", "--secondaryPeaksCount", dest="secondaryPeaksCount", type=int, default=10,
                            help="Number of peaks found in the second cross-correlation run that are selected for alignment creation "
                                 "(for each query molecule and its peak from the first cross-correlation). ")

        parser.add_argument("-ma", "--secondaryMargin", dest="secondaryMargin", type=int, default=100000,
                            help="The number of base pairs by which the peak from initial cross-correlation "
                                 "seeding is extended in both directions to serve as an input "
                                 "for the second cross-correlation run.")

        parser.add_argument("-pt", "--peakHeightThreshold", dest="peakHeightThreshold", type=float, default=10,
                            help="Minimum second cross-correlation peak height to qualify for aligned pairs search.")

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

        parser.add_argument("-ms", "--minScore", dest="minScore", type=int, default=20000,
                            help="Minimum score of a segment/alignment.")

        parser.add_argument("-mss", "--minSubsequentScore", dest="minSubsequentScore", type=int, default=None,
                            help="Minimum alignment score in subsequent passes.")

        parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=1200,
                            help="Alignment segments can be split into two if their score drops below this threshold.")

        parser.add_argument("-diff", "--maxDifference", dest="maxDifference", type=int, default=100000,
                            help="Multiple alignments of the same query will be joined if difference between their "
                                 "reference positions is less or equal this parameter.")

        parser.add_argument("-D", "--diagnostics", dest="diagnosticsEnabled", action="store_true",
                            help="Draws cross-correlation and alignment plots. When used, 'outputFile' parameter "
                                 "is required. When 'benchmarkAlignmentFile' is provided, alignment plots will allow "
                                 "to compare both alignments if they are overlapping.")

        parser.add_argument("-a", "--benchmarkAlignment", dest="benchmarkAlignmentFile", type=argparse.FileType("r"),
                            default=None,
                            help="XMAP file containing alignments from other source, to be used with 'diagnostics' "
                                 "option.")

        parser.add_argument("-c", "--cpus", dest="numberOfCpus", type=int, default=None,
                            help="Number of CPUs to use. The default is all available CPUs.")

        parser.add_argument("-pb", "--disableProgressBar", dest="disableProgressBar", action="store_true",
                            help="Disables the progress bar.")

        parser.add_argument("-er", "--endReachingScore", dest="endReachingScore", type=int, default=10000,
                            help="Score value given to an alignment reaching query start or end.")

        parser.add_argument("-sc", "--segmentCombinePenalty", dest="segmentCombinePenalty", type=int, default=5000,
                            help="Constant component of penalty for combining segment.")

        parser.add_argument("-sj", "--segmentJoinMultiplier", dest="segmentJoinMultiplier", type=float, default=0.1,
                            help="Multiplier applied to segment sequentiality scores.")

        parser.add_argument("-ss", "--sequentialityScore", dest="sequentialityScore", type=int, default=2,
                            help="Segment sequentiality scoring function version.")
        
        args = parser.parse_args(args)
        return args  # type: ignore
