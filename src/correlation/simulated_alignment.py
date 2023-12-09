from typing import List

from src.correlation.bionano_alignment import BionanoAlignment
from src.diagnostic.benchmark_alignment import BenchmarkAlignment, BenchmarkAlignedPair


class SimulatedAlignment(BenchmarkAlignment):
    def __init__(self, alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, reverseStrand, confidence,
                 cigarString, queryLength, referenceLength, alignedPairs: List[BenchmarkAlignedPair]) -> None:
        self.alignmentId = alignmentId
        self.queryId = queryId
        self.referenceId = refId
        self.queryStartPosition = queryStart
        self.queryEndPosition = queryEnd
        self.referenceStartPosition = refStart
        self.referenceEndPosition = refEnd
        self.reverseStrand = reverseStrand
        self.confidence = confidence
        self.cigarString = cigarString
        self.queryLength = queryLength
        self.referenceLength = referenceLength
        self.alignedPairs = alignedPairs

    @staticmethod
    def parse(alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, reverseStrand, confidence,
              cigarString, queryLength, referenceLength, alignment: List[BenchmarkAlignedPair]):
        return BionanoAlignment(
            alignmentId,
            int(queryId),
            int(refId),
            int(queryStart),
            int(queryEnd),
            int(refStart),
            int(refEnd),
            reverseStrand,
            confidence,
            cigarString,
            int(queryLength),
            int(referenceLength),
            alignment)
