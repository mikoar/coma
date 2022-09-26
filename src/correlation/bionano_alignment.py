from typing import List

from src.diagnostic.xmap_alignment import XmapAlignment, XmapAlignedPair


class BionanoAlignment(XmapAlignment):
    def __init__(self, alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, reverseStrand, confidence,
                 cigarString, queryLength, referenceLength, alignedPairs: List[XmapAlignedPair]) -> None:
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
              cigarString, queryLength, referenceLength, alignment: List[XmapAlignedPair]):
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

    @property
    def expectedQueryMoleculeStart(self):
        return self.referenceStartPosition - self.__queryStartPositionDisregardingOrientation

    @property
    def expectedQueryMoleculeEnd(self):
        return self.referenceEndPosition + self.queryLength - self.__queryEndPositionDisregardingOrientation

    @property
    def queryReferenceAlignmentLengthDifference(self):
        """Alignment length on reference and query sequences may differ due to indels and molecule stretch"""
        return (self.queryAlignmentLength()) - (self.referenceAlignmentLength())

    def referenceAlignmentLength(self):
        return abs(self.referenceEndPosition - self.referenceStartPosition)

    def queryAlignmentLength(self):
        return abs(self.__queryEndPositionDisregardingOrientation - self.__queryStartPositionDisregardingOrientation)

    @property
    def __queryStartPositionDisregardingOrientation(self):
        return self.queryLength - self.queryStartPosition if self.reverseStrand else self.queryStartPosition

    @property
    def __queryEndPositionDisregardingOrientation(self):
        return self.queryLength - self.queryEndPosition if self.reverseStrand else self.queryEndPosition
