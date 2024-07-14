from __future__ import annotations

import math
from typing import Iterable, List

from src.alignment.segments import AlignmentSegment


class SegmentChainer:
    def __init__(self, sequentialityScorer: SequentialityScorer):
        self.sequentialityScorer = sequentialityScorer

    def chain(self, segments: Iterable[AlignmentSegment]):
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
        for i, currentSegment in enumerate(preOrderedNonEmptySegments):
            cumulatedScore[i] = 0
            for j, previousSegment in enumerate(preOrderedNonEmptySegments[:i]):
                currentScore = cumulatedScore[j] + self.sequentialityScorer.getScore(
                    previousSegment, currentSegment)
                if currentScore > cumulatedScore[i]:
                    cumulatedScore[i] = currentScore
                    previousSegmentIndexes[i] = j
            cumulatedScore[i] += currentSegment.segmentScore
            if cumulatedScore[i] > cumulatedScore[bestPreviousSegmentIndex]:
                bestPreviousSegmentIndex = i
        result = [preOrderedNonEmptySegments[bestPreviousSegmentIndex]]
        while (bestPreviousSegmentIndex := previousSegmentIndexes[bestPreviousSegmentIndex]) is not None:
            result.insert(0, preOrderedNonEmptySegments[bestPreviousSegmentIndex])
        return result + emptySegments


class SequentialityScorer:
    def __init__(self, segmentJoinMultiplier: float, sequentialityScore: int):
        self.segmentJoinMultiplier = segmentJoinMultiplier
        self.sequentialityScore = sequentialityScore

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

        def calcScore(referDist, queryDist, seqScore):
            distSum = referDist + queryDist
            absDistSum = abs(referDist) + abs(queryDist)
            distDiff = referDist - queryDist
            if seqScore == 0:
                return (distSum ** 2 + distDiff ** 2) / max(abs(distSum), abs(distDiff), 1)
            else:
                return (absDistSum ** 2 + distDiff ** 2) / max(absDistSum + abs(distDiff), 1)

        return - self.segmentJoinMultiplier * calcScore(referenceDistance, queryDistance, self.sequentialityScore)
