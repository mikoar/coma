from __future__ import annotations

import argparse
import sys
from typing import NamedTuple, TextIO, List


class Args(NamedTuple):
    referenceFile: TextIO
    queryFile: TextIO
    benchmarkAlignmentFile: TextIO
    outputFile: TextIO
    primaryResolution: int
    primaryBlur: int
    secondaryResolution: int
    secondaryBlur: int
    adjustment: int
    referenceIds: List[int]
    queryIds: List[int]
    numberOfCpus: int | None
    maxDistance: int
    peakHeightThreshold: float
    perfectMatchScore: int
    scoreMultiplier: int
    unmatchedPenalty: int
    minScore: int
    breakSegmentThreshold: int
    diagnosticsEnabled: bool

    @staticmethod
    def parse() -> Args:
        parser = argparse.ArgumentParser(description="Optical map aligner.")
        parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"), required=True)
        parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"), required=True)
        parser.add_argument("-a", "--benchmarkAlignment", dest="benchmarkAlignmentFile", type=argparse.FileType("r"),
                            default=None)
        parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"),
                            default=sys.stdout)
        parser.add_argument("-r1", "--primaryResolution", dest="primaryResolution", type=int, default=460)
        parser.add_argument("-b1", "--primaryBlur", dest="primaryBlur", type=int, default=2)
        parser.add_argument("-r2", "--secondaryResolution", dest="secondaryResolution", type=int, default=200)
        parser.add_argument("-b2", "--secondaryBlur", dest="secondaryBlur", type=int, default=2)
        parser.add_argument("-ma", "--adjustment", dest="adjustment", type=int, default=8000)
        parser.add_argument("-rId", "--referenceIDs", dest="referenceIds", type=int, nargs="*")
        parser.add_argument("-qId", "--queryIDs", dest="queryIds", type=int, nargs="*")
        parser.add_argument("-c", "--cpus", dest="numberOfCpus", type=int, default=None)
        parser.add_argument("-d", "--maxDistance", dest="maxDistance", type=int, default=1000)
        parser.add_argument("-pt", "--peakHeightThreshold", dest="peakHeightThreshold", type=float, default=15)
        parser.add_argument("-sp", "--perfectMatchScore", dest="perfectMatchScore", type=int, default=800)
        parser.add_argument("-sm", "--scoreMultiplier", dest="scoreMultiplier", type=float, default=1.)
        parser.add_argument("-su", "--unmatchedPenalty", dest="unmatchedPenalty", type=int, default=-100)
        parser.add_argument("-ms", "--minScore", dest="minScore", type=int, default=1600)
        parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=600)
        parser.add_argument("-D", "--diagnostics", dest="diagnosticsEnabled", action="store_true")
        args = parser.parse_args()
        return args  # type: ignore
