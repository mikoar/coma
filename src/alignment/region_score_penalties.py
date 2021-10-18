from __future__ import annotations

from src.alignment.aligned_pair import AlignedPair


class RegionScorePenalties:
    def __init__(self, unmatchedLabelPenalty: int = 5000, distancePenaltyExponent: float = 1.2) -> None:
        self.unmatchedLabelPenalty = unmatchedLabelPenalty
        self.distancePenaltyExponent = distancePenaltyExponent

    def getUnmatchedLabelPenalty(self, previousPair: AlignedPair | None, currentPair: AlignedPair):
        if not previousPair:
            return 0
        else:
            return self.unmatchedLabelPenalty * (currentPair.getFalseNegativesCount(previousPair)
                                                 + currentPair.getFalsePositivesCount(previousPair))

    def getDistancePenalty(self, previousPair: AlignedPair | None, currentPair: AlignedPair):
        if not previousPair:
            return 0
        return abs(previousPair.queryShift - currentPair.queryShift) ** self.distancePenaltyExponent
