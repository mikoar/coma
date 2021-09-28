from __future__ import annotations
from typing import Iterable, List, NamedTuple, Union

from .aligner import AlignmentResult, AlignedPair


class AlignmentSegment(NamedTuple):
    start: int
    end: int
    score: float


class AlignmentScorer:
    def getSegmentWithMaxScore(self, alignmentScores: List[float]):
        start = end = 0
        currentScore = 0
        segmentWithMaxScore = AlignmentSegment(0, 0, 0)

        alignmentEnd = len(alignmentScores) - 1
        while end <= alignmentEnd:
            if currentScore + alignmentScores[end] > 0:
                currentScore += alignmentScores[end]
                if currentScore > segmentWithMaxScore.score:
                    segmentWithMaxScore = AlignmentSegment(start, end, currentScore)
                end += 1
            else:
                start = end = end + 1
                currentScore = 0

        return segmentWithMaxScore


class RegionScorer:
    def __init__(self, penalties: RegionScorePenalties, perfectMatchScore: int = 10000) -> None:
        self.perfectMatchScore = perfectMatchScore
        self.penalties = penalties

    def getRegionScores(self, alignmentResult: AlignmentResult) -> Iterable[Union[int, float]]:
        previousPair: AlignedPair | None = None
        for pair in alignmentResult.alignedPairs:
            yield (self.perfectMatchScore
                   - self.penalties.getUnmatchedLabelPenalty(previousPair, pair)
                   - self.penalties.getDistancePenalty(previousPair, pair))
            previousPair = pair


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
        return 0.
