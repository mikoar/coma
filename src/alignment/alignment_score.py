from __future__ import annotations

from typing import List, NamedTuple

from src.alignment.aligned_pair import AlignedPair


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
