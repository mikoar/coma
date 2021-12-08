from typing import NamedTuple, List


class RefAlignedPair(NamedTuple):
    referencePositionIndex: int
    queryPositionIndex: int


class BionanoAlignment:
    def __init__(self, alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, reverseStrand, confidence,
                 hitEnum, queryLength, referenceLength, alignedPairs: List[RefAlignedPair]) -> None:
        self.alignmentId = alignmentId
        self.queryId = queryId
        self.referenceId = refId
        self.__queryStartPositionRelativeToStrand = queryStart
        self.__queryEndPositionRelativeToStrand = queryEnd
        self.referenceStartPosition = refStart
        self.referenceEndPosition = refEnd
        self.reverseStrand = reverseStrand
        self.confidence = confidence
        self.cigarString = hitEnum
        self.queryLength = queryLength
        self.referenceLength = referenceLength
        self.alignedPairs = alignedPairs

    @staticmethod
    def parse(alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, orientation, confidence,
              hitEnum, queryLength, referenceLength, alignment: str):
        return BionanoAlignment(
            alignmentId,
            int(queryId),
            int(refId),
            int(queryStart),
            int(queryEnd),
            int(refStart),
            int(refEnd),
            orientation == "-",
            confidence,
            hitEnum,
            int(queryLength),
            int(referenceLength),
            list(map(lambda pair: RefAlignedPair(*pair.split(',')), alignment[:-1].replace('(', '').split(')'))))

    @property
    def expectedPeakPosition(self):
        expectedQueryMoleculeStart = self.expectedQueryMoleculeStart
        moleculeLength = self.expectedQueryMoleculeEnd - expectedQueryMoleculeStart
        return expectedQueryMoleculeStart + moleculeLength / 2

    @property
    def queryStartPosition(self):
        return self.queryLength - self.__queryStartPositionRelativeToStrand if self.reverseStrand else self.__queryStartPositionRelativeToStrand

    @property
    def queryEndPosition(self):
        return self.queryLength - self.__queryEndPositionRelativeToStrand if self.reverseStrand else self.__queryEndPositionRelativeToStrand

    @property
    def expectedQueryMoleculeStart(self):
        return self.referenceStartPosition - self.queryStartPosition - self.queryReferenceAlignmentLengthDifference

    @property
    def expectedQueryMoleculeEnd(self):
        return self.referenceEndPosition + self.queryLength - self.queryEndPosition + self.queryReferenceAlignmentLengthDifference

    @property
    def queryReferenceAlignmentLengthDifference(self):
        """Alignment length on reference and query sequences may differ due to insertions and deletions"""
        return (self.queryAlignmentLength()) - (self.referenceAlignmentLength())

    def referenceAlignmentLength(self):
        return abs(self.referenceEndPosition - self.referenceStartPosition)

    def queryAlignmentLength(self):
        return abs(self.queryEndPosition - self.queryStartPosition)
