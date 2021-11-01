from typing import NamedTuple, List


class RefAlignedPair(NamedTuple):
    referencePositionIndex: int
    queryPositionIndex: int


class Alignment:
    def __init__(self, alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, reverseStrand, confidence,
                 queryLength, alignedPairs: List[RefAlignedPair]) -> None:
        self.alignmentId = alignmentId
        self.queryId = queryId
        self.referenceId = refId
        self.__queryStartPositionRelativeToStrand = queryStart
        self.__queryEndPositionRelativeToStrand = queryEnd
        self.referenceAlignmentStartPosition = refStart
        self.referenceAlignmentEndPosition = refEnd
        self.reverseStrand = reverseStrand
        self.confidence = confidence
        self.queryLength = queryLength
        self.alignedPairs = alignedPairs

    @staticmethod
    def parse(alignmentId, queryId, refId, queryStart, queryEnd, refStart, refEnd, orientation, confidence,
              queryLength, alignment: str):
        return Alignment(
            alignmentId,
            int(queryId),
            int(refId),
            int(queryStart),
            int(queryEnd),
            int(refStart),
            int(refEnd),
            orientation == "-",
            confidence,
            int(queryLength),
            list(map(lambda pair: RefAlignedPair(*pair.split(',')), alignment[:-1].replace('(', '').split(')'))))

    @property
    def expectedPeakPosition(self):
        expectedQueryMoleculeStart = self.expectedQueryMoleculeStart
        moleculeLength = self.expectedQueryMoleculeEnd - expectedQueryMoleculeStart
        return expectedQueryMoleculeStart + moleculeLength / 2

    @property
    def queryAlignmentStartPosition(self):
        return self.queryLength - self.__queryStartPositionRelativeToStrand if self.reverseStrand else self.__queryStartPositionRelativeToStrand

    @property
    def queryAlignmentEndPosition(self):
        return self.queryLength - self.__queryEndPositionRelativeToStrand if self.reverseStrand else self.__queryEndPositionRelativeToStrand

    @property
    def expectedQueryMoleculeStart(self):
        return self.referenceAlignmentStartPosition - self.queryAlignmentStartPosition - self.queryReferenceAlignmentLengthDifference

    @property
    def expectedQueryMoleculeEnd(self):
        return self.referenceAlignmentEndPosition + self.queryLength - self.queryAlignmentEndPosition + self.queryReferenceAlignmentLengthDifference

    @property
    def queryReferenceAlignmentLengthDifference(self):
        """Alignment length on reference and query sequences may differ due to insertions and deletions"""
        return (self.queryAlignmentLength()) - (self.referenceAlignmentLength())

    def referenceAlignmentLength(self):
        return abs(self.referenceAlignmentEndPosition - self.referenceAlignmentStartPosition)

    def queryAlignmentLength(self):
        return abs(self.queryAlignmentEndPosition - self.queryAlignmentStartPosition)
