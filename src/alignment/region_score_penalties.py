from __future__ import annotations

from abc import ABC, abstractmethod

from src.alignment.aligned_pair import AlignedPair


class RegionScorePenalty(ABC):
    @abstractmethod
    def getPenalty(self, previousPair: AlignedPair | None, currentPair: AlignedPair) -> float:
        pass


class UnmatchedLabelPenalty(RegionScorePenalty):
    def __init__(self, penalty: int):
        self.penalty = penalty

    def getPenalty(self, previousPair: AlignedPair | None, currentPair: AlignedPair) -> float:
        if not previousPair:
            return 0
        else:
            return self.penalty * (currentPair.getFalseNegativesCount(previousPair) +
                                   currentPair.getFalsePositivesCount(previousPair))


class DistancePenalty(RegionScorePenalty):
    def __init__(self, multiplier: int = 2):
        self.multiplier = multiplier

    def getPenalty(self, previousPair: AlignedPair | None, currentPair: AlignedPair) -> float:
        if not previousPair:
            return 0
        return abs(previousPair.queryShift - currentPair.queryShift) * self.multiplier
