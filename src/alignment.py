class Alignment:
    def __init__(self, id, queryId, refId, queryStart, queryEnd, refStart, refEnd, orientation, confidence, queryLength) -> None:
        self.id = id
        self.queryId = int(queryId)
        self.referenceId = int(refId)
        self.__queryStartPositionRelativeToStrand = int(queryStart)
        self.__queryEndPositionRelativeToStrand = int(queryEnd)
        self.referenceAlignmentStartPosition = int(refStart)
        self.referenceAlignmentEndPosition = int(refEnd)
        self.reverseStrand = orientation == "-"
        self.confidence = confidence
        self.queryLength = int(queryLength)

    @property
    def expectedPeakPosition(self):
        expectedQueryStart = self.queryAlignmentStartPosition
        length = self.queryAlignmentEndPosition - expectedQueryStart
        return expectedQueryStart + length / 2

    @property
    def queryAlignmentStartPosition(self):
        return self.queryLength - self.__queryStartPositionRelativeToStrand if self.reverseStrand else self.__queryStartPositionRelativeToStrand

    @property
    def queryAlignmentEndPosition(self):
        return self.queryLength - self.__queryEndPositionRelativeToStrand if self.reverseStrand else self.__queryEndPositionRelativeToStrand

    @property
    def expectedQueryMoleculeStart(self):
        return self.referenceAlignmentStartPosition - self.queryAlignmentStartPosition - self.__queryReferenceAlignmentLengthDifference

    @property
    def expectedQueryMoleculeEnd(self):
        return self.referenceAlignmentEndPosition + self.queryLength - self.queryAlignmentEndPosition + self.__queryReferenceAlignmentLengthDifference

    @property
    def __queryReferenceAlignmentLengthDifference(self):
        """Alignment length on reference and query sequences may differ due to insertions and deletions"""
        return (self.queryAlignmentEndPosition - self.queryAlignmentStartPosition) - (self.referenceAlignmentEndPosition - self.referenceAlignmentStartPosition)
