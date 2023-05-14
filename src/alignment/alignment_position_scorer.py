from typing import List

from src.alignment.alignment_position import AlignmentPosition


class AlignmentPositionScorer:
    def __init__(self, perfectMatchScore: int,
                 distancePenaltyMultiplier: float,
                 unmatchedPenalty: int):
        self.perfectMatchScore = perfectMatchScore
        self.distancePenaltyMultiplier = distancePenaltyMultiplier
        self.unmatchedPenalty = unmatchedPenalty

    def getScoredPositions(self, positions: List[AlignmentPosition]):
        return [p.getScoredPosition(self.perfectMatchScore, self.distancePenaltyMultiplier, self.unmatchedPenalty) for p
                in positions]
