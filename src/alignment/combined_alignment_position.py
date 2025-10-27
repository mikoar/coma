from __future__ import annotations

from abc import ABC

from src.alignment.alignment_position import AlignmentPosition, AlignedPair, NotAlignedReferencePosition, NotAlignedQueryPosition
# from src.alignment_position_scorer import AlignmentPositionScorer


class AlignmentPositionScorer:
    def __init__(self, perfectMatchScore: int,
                 distancePenaltyMultiplier: float,
                 unmatchedPenalty: int):
        self.perfectMatchScore = perfectMatchScore
        self.distancePenaltyMultiplier = distancePenaltyMultiplier
        self.unmatchedPenalty = unmatchedPenalty

class StartEndPosition(AlignmentPosition):
    def __init__(self, queryPosition: int):
        self.queryPosition = queryPosition

class CombinedAlignmentPosition(ABC):
    score: float
    queryPosition: int
    referenceStart: int
    alignmentPosition: AlignmentPosition
    def __init__(self, alignmentPosition: AlignmentPosition, score: float, referenceStart: int):
        self.alignmentPosition = alignmentPosition
        self.score = score
        self.referenceStart = referenceStart
        self.setQueryPosition(alignmentPosition)

    def __lt__(self, other: CombinedAlignmentPosition):
        return self.queryPosition < other.queryPosition
    
    def lessOnBothMaps(self, other: CombinedAlignmentPosition):
        return self.queryPosition < other.queryPosition and self.referencePosition() < other.referencePosition()
    
    def referencePosition(self):
        return self.referenceStart + self.queryPosition


class CombinedStartEndPosition(CombinedAlignmentPosition):
    def setQueryPosition(self, position: StartEndPosition):
        self.queryPosition = position.queryPosition

class CombinedAlignedPair(CombinedAlignmentPosition):
    def setQueryPosition(self, pair: AlignedPair):
        self.queryPosition = pair.query.position

class CombinedNAQueryPosition(CombinedAlignmentPosition):
    def setQueryPosition(self, queryPosition: NotAlignedQueryPosition):
        self.queryPosition = queryPosition.query.position

class CombinedNAReferPosition(CombinedAlignmentPosition):
    def setQueryPosition(self, referencePosition: NotAlignedReferencePosition):
        self.queryPosition = referencePosition.reference.position - self.referenceStart
    




