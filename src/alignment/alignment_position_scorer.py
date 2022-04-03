from typing import List

from src.alignment.alignment_position import AlignmentPosition


class AlignmentPositionScorer:
    def __init__(self, perfectMatchScore: int,
                 scoreMultiplier: int,
                 unmatchedPenalty: int):
        self.perfectMatchScore = perfectMatchScore
        self.scoreMultiplier = scoreMultiplier
        self.unmatchedPenalty = unmatchedPenalty

    def getScoredPositions(self, positions: List[AlignmentPosition]):
        return [p.getScoredPosition(self.perfectMatchScore, self.scoreMultiplier, self.unmatchedPenalty) for p in
                positions]
