from __future__ import annotations

import math
from typing import Iterable, List

from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import PositionWithSiteId


class SegmentChainer:
    def __init__(self, sequentialityScorer: SequentialityScorer):
        self.sequentialityScorer = sequentialityScorer

    def chain(self, segments: Iterable[AlignmentSegment], queryPositions: list[PositionWithSiteId]):
        def initialOrderingKey(segment: AlignmentSegment):
            return segment.startPosition.reference.position + segment.endPosition.reference.position \
                   + segment.startPosition.query.position + segment.endPosition.query.position

        emptySegments = [s for s in segments if s.empty]
        preOrderedNonEmptySegments = sorted((s for s in segments if not s.empty), key=initialOrderingKey)
        if not preOrderedNonEmptySegments:
            return emptySegments

        cumulatedScore = [-math.inf] * len(preOrderedNonEmptySegments)
        previousSegmentIndexes: List[int | None] = [None] * len(preOrderedNonEmptySegments)
        bestPreviousSegmentIndex = 0
        bestPreviousTotalScore = -math.inf
        for i, currentSegment in enumerate(preOrderedNonEmptySegments):
            cumulatedScore[i] = self.sequentialityScorer.getStartScore(queryPositions[0], currentSegment)
            for j, previousSegment in enumerate(preOrderedNonEmptySegments[:i]):
                currentScore = cumulatedScore[j] + self.sequentialityScorer.getScore(
                    previousSegment, currentSegment)
                if currentScore > cumulatedScore[i]:
                    cumulatedScore[i] = currentScore
                    previousSegmentIndexes[i] = j
            cumulatedScore[i] += currentSegment.segmentScore
            currentTotalScore = cumulatedScore[i] + self.sequentialityScorer.getEndScore(queryPositions[-1], currentSegment)
            if currentTotalScore > bestPreviousTotalScore:
                bestPreviousSegmentIndex = i
                bestPreviousTotalScore = currentTotalScore
        result = [preOrderedNonEmptySegments[bestPreviousSegmentIndex]]
        while (bestPreviousSegmentIndex := previousSegmentIndexes[bestPreviousSegmentIndex]) is not None:
            result.insert(0, preOrderedNonEmptySegments[bestPreviousSegmentIndex])
        return result + emptySegments


class SequentialityScorer:
    def __init__(self, segmentJoinMultiplier: float, sequentialityScore: int, segmentBreakPenalty: int=0, secondaryMargin: int=100000):
        self.segmentJoinMultiplier = segmentJoinMultiplier
        self.sequentialityScore = sequentialityScore
        self.segmentBreakPenalty = segmentBreakPenalty
        self.secondaryMargin = secondaryMargin

    def calcScore(self, referDist, queryDist, seqScore):
        distSum = referDist + queryDist
        absDistSum = abs(referDist) + abs(queryDist)
        distDiff = referDist - queryDist
        if seqScore == 0:
            return (distSum ** 2 + distDiff ** 2) / max(abs(distSum), abs(distDiff), 1)
        else:
            return (absDistSum ** 2 + distDiff ** 2) / max(absDistSum + abs(distDiff), 1)

    def getScore(self, previousSegment: AlignmentSegment, currentSegment: AlignmentSegment):
        queryLength = min(abs(currentSegment.endPosition.query.position - currentSegment.startPosition.query.position),
                          abs(previousSegment.endPosition.query.position - previousSegment.startPosition.query.position))
        referenceDistance = currentSegment.startPosition.reference.position - previousSegment.endPosition.reference.position
        referenceLength = min(
            currentSegment.endPosition.reference.position - currentSegment.startPosition.reference.position,
            previousSegment.endPosition.reference.position - previousSegment.startPosition.reference.position)

        queryDistance = previousSegment.endPosition.query.position - currentSegment.startPosition.query.position \
            if currentSegment.reverse \
            else currentSegment.startPosition.query.position - previousSegment.endPosition.query.position

        if min(referenceLength + 2 * referenceDistance, queryLength + 2 * queryDistance) < 0:
            return -math.inf

        
        distScore = self.calcScore(referenceDistance, queryDistance, self.sequentialityScore)
        return - self.segmentBreakPenalty - self.segmentJoinMultiplier * distScore

    def getStartScore(self, queryStart: PositionWithSiteId, currentSegment: AlignmentSegment):
        queryDistance = min(self.secondaryMargin, currentSegment.startPosition.query.position - queryStart.position)
        assert queryDistance>=0
        distScore = self.calcScore(0, queryDistance, self.sequentialityScore)
        return - self.segmentBreakPenalty - self.segmentJoinMultiplier * distScore

    def getEndScore(self, queryEnd: PositionWithSiteId, currentSegment: AlignmentSegment):
        queryDistance = min(self.secondaryMargin, queryEnd.position - currentSegment.endPosition.query.position)
        assert queryDistance>=0
        distScore = self.calcScore(0, queryDistance, self.sequentialityScore)
        return - self.segmentBreakPenalty - self.segmentJoinMultiplier * distScore


