from __future__ import annotations

# from abc import ABC

from src.alignment.alignment_position import AlignmentPosition, AlignedPair, NotAlignedReferencePosition, NotAlignedQueryPosition
# from src.alignment_position_scorer import AlignmentPositionScorer



class MockPosition(AlignmentPosition):
    def __init__(self, queryPosition: int, start: bool):
        self.queryPosition = queryPosition
        self.start = start
    def getScoredPosition(self, perfectMatchScore: int, distancePenaltyMultiplier: float,
                          unmatchedPenalty: int):
        pass
    def absolutePosition(self):
        pass

class ChainedAlignmentPosition:
    def __init__(self, alignmentPosition: AlignmentPosition, score: float, referenceStart: int):
        self.alignmentPosition = alignmentPosition
        self.score = score
        self.queryShift = referenceStart
        if isinstance(alignmentPosition, MockPosition):
            self.queryPosition = alignmentPosition.queryPosition
            self.referenceIdx = -1
            self.queryIdx = -1
        elif isinstance(alignmentPosition, NotAlignedReferencePosition):
            self.queryPosition = alignmentPosition.reference.position - self.queryShift
            self.referenceIdx = alignmentPosition.reference.siteId
            self.queryIdx = -1
        elif isinstance(alignmentPosition, NotAlignedQueryPosition):
            self.queryPosition = alignmentPosition.query.position
            self.referenceIdx = -1
            self.queryIdx = alignmentPosition.query.siteId
        else:
            self.queryPosition = alignmentPosition.query.position
            self.referenceIdx = alignmentPosition.reference.siteId
            self.queryIdx = alignmentPosition.query.siteId

    def __lt__(self, other: ChainedAlignmentPosition):
        return self.queryPosition < other.queryPosition
    
    def lessOnBothMaps(self, other: ChainedAlignmentPosition):
        return self.queryPosition < other.queryPosition and self.referencePosition() < other.referencePosition()
    
    def referencePosition(self):
        return self.queryShift + self.queryPosition
    
    def queryPosition(self):
        return self.queryPosition

    def isNARefer(self):
        return isinstance(self.alignmentPosition, NotAlignedReferencePosition)

    def isNAQuery(self):
        return isinstance(self.alignmentPosition, NotAlignedQueryPosition)

    def isAlPair(self):
        return isinstance(self.alignmentPosition, AlignedPair)

    def isStart(self):
        return isinstance(self.alignmentPosition, MockPosition) and self.alignmentPosition.start

    def isEnd(self):
        return isinstance(self.alignmentPosition, MockPosition) and not self.alignmentPosition.start





