from __future__ import annotations

import argparse
import sys
from typing import NamedTuple, TextIO, List


class Args(NamedTuple):
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO
    primaryResolution: int
    primaryBlur: int
    secondaryResolution: int
    secondaryBlur: int
    secondaryMargin: int
    referenceIds: List[int]
    queryIds: List[int]
    numberOfCpus: int | None
    minPeakDistance: int
    maxDistance: int
    peakHeightThreshold: float
    perfectMatchScore: int
    distancePenaltyMultiplier: int
    unmatchedPenalty: int
    minScore: int
    breakSegmentThreshold: int
    diagnosticsEnabled: bool
    benchmarkAlignmentFile: TextIO

    @staticmethod
    def parse() -> Args:
        parser = argparse.ArgumentParser(description="Optical map aligner.")

        parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"), required=True,
                            help="Reference optical map in CMAP format file path.")

        parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"), required=True,
                            help="Query optical map in CMAP format file path.")

        parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"),
                            default=sys.stdout,
                            help="XMAP output file path. Stdout is used if omitted.")

        parser.add_argument("-r1", "--primaryResolution", dest="primaryResolution", type=int, default=460,
                            help="Scaling factor used to reduce the size of the vectorized form of the optical map "
                                 "in the initial cross-correlation seeding step.")

        parser.add_argument("-b1", "--primaryBlur", dest="primaryBlur", type=int, default=2,
                            help="Extends each label in the vectorized form of the optical map in both directions "
                                 "by given number of positions in the initial cross-correlation seeding step "
                                 "in order to increase the chance of overlap. "
                                 "Final width of each label is equal to 2b + 1.")

        parser.add_argument("-r2", "--secondaryResolution", dest="secondaryResolution", type=int, default=200,
                            help="Scaling factor used to reduce the size of the vectorized form of the optical map "
                                 "in the second cross-correlation run.")

        parser.add_argument("-b2", "--secondaryBlur", dest="secondaryBlur", type=int, default=2,
                            help="Extends each label in the vectorized form of the optical map in both directions "
                                 "by given number of positions in the second cross-correlation run "
                                 "in order to increase the chance of overlap. "
                                 "Final width of each label is equal to 2b + 1.")

        parser.add_argument("-ma", "--secondaryMargin", dest="secondaryMargin", type=int, default=8000,
                            help="The number of base pairs by which the peak from initial cross-correlation "
                                 "seeding is extended in both directions to serve as an input "
                                 "for the second cross-correlation run.")

        parser.add_argument("-rId", "--referenceIDs", dest="referenceIds", type=int, nargs="*",
                            help="CMapId(s) of reference molecules to be used. Takes all if omitted.")

        parser.add_argument("-qId", "--queryIDs", dest="queryIds", type=int, nargs="*",
                            help="CMapId(s) of query molecules to be used. Takes all if omitted.")

        parser.add_argument("-c", "--cpus", dest="numberOfCpus", type=int, default=None,
                            help="Number of CPUs to use. The default is all CPUs.")

        parser.add_argument("-md", "--minPeakDistance", dest="minPeakDistance", type=int, default=20000,
                            help="Minimum distance between peaks identified in the initial cross-correlation. "
                                 "For more details see parameter distance of scipy.signal._peak_finding.find_peaks.")

        parser.add_argument("-d", "--maxDistance", dest="maxDistance", type=int, default=1000,
                            help="Maximum distance between aligned pairs relatively to the cross-correlation lag.")

        parser.add_argument("-pt", "--peakHeightThreshold", dest="peakHeightThreshold", type=float, default=15,
                            help="Minimum second cross-correlation peak height to qualify for aligned pairs search.")

        parser.add_argument("-sp", "--perfectMatchScore", dest="perfectMatchScore", type=int, default=800,
                            help="Score value given to an aligned pair with 0 distance between reference and query "
                                 "positions.")

        parser.add_argument("-dp", "--distancePenaltyMultiplier", dest="distancePenaltyMultiplier", type=float,
                            default=1., help="Multiplier applied to the distance between reference and query positions "
                                             "of an aligned pair that reduces the pair's score.")

        parser.add_argument("-su", "--unmatchedPenalty", dest="unmatchedPenalty", type=int, default=-100,
                            help="Penalty to a segment score for each unpaired reference or query position.")

        parser.add_argument("-ms", "--minScore", dest="minScore", type=int, default=1000,
                            help="Minimum score of a segment.")

        parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=1200,
                            help="Alignment segments can be split into two if their score drops below this threshold.")

        parser.add_argument("-D", "--diagnostics", dest="diagnosticsEnabled", action="store_true",
                            help="Draws cross-correlation and alignment plots. When used, 'outputFile' parameter "
                                 "is required. When 'benchmarkAlignmentFile' is provided, alignment plots will allow "
                                 "to compare both alignments if they are overlapping.")

        parser.add_argument("-a", "--benchmarkAlignment", dest="benchmarkAlignmentFile", type=argparse.FileType("r"),
                            default=None,
                            help="XMAP file containing alignments from other source, to be used with 'diagnostics' "
                                 "option.")

        args = parser.parse_args()
        return args  # type: ignore
