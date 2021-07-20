class Alignment:
    def __init__(self, id, queryId, refId, queryStart, queryEnd, refStart, refEnd, orientation, confidence, queryLength) -> None:
        self.id = id
        self.queryId = int(queryId)
        self.referenceId = int(refId)
        self.queryStartPosition = int(queryStart)
        self.queryEndPosition = int(queryEnd)
        self.refStartPosition = int(refStart)
        self.refEndPosition = int(refEnd)
        self.reverseStrand = orientation == "-"
        self.confidence = confidence
        self.queryLength = int(queryLength)

    @property
    def expectedPeakPosition(self):
        expectedQueryStart = self.expectedQueryStart
        length = self.expectedQueryEnd - expectedQueryStart
        return expectedQueryStart + length / 2

    @property
    def expectedQueryStart(self):
        return self.refStartPosition - self.queryStartPosition - self.__queryReferenceAlignmentLengthDifference

    @property
    def expectedQueryEnd(self):
        return self.refEndPosition + self.queryLength - self.queryEndPosition + self.__queryReferenceAlignmentLengthDifference

    @property
    def __queryReferenceAlignmentLengthDifference(self):
        """Alignment length on reference and query sequences may differ due to insertions and deletions"""
        return (self.queryEndPosition - self.queryStartPosition) - (self.refEndPosition - self.refStartPosition)
