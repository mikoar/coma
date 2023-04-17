from __future__ import annotations

import math
from typing import Iterable, List

from src.alignment.segments import AlignmentSegment


class SegmentChainer:
    def chain(self, segments: Iterable[AlignmentSegment]):
        def initialOrderingKey(segment: AlignmentSegment):
            return segment.startPosition.reference.position + segment.endPosition.reference.position \
                   + segment.startPosition.query.position + segment.endPosition.query.position

        segments = sorted(segments, key=initialOrderingKey)
        cumulatedScore = [-math.inf] * len(segments)
        previousSegmentIndexes: List[int | None] = [None] * len(segments)
        bestPreviousSegmentIndex = 0
        for i, currentSegment in enumerate(segments):
            cumulatedScore[i] = 0
            for j, previousSegment in enumerate(segments[:i]):
                currentScore = cumulatedScore[j] + self.__consecutivenessScore(previousSegment,
                                                                               currentSegment)
                if currentScore > cumulatedScore[i]:
                    cumulatedScore[i] = currentScore
                    previousSegmentIndexes[i] = j
            cumulatedScore[i] += currentSegment.segmentScore
            if cumulatedScore[i] > cumulatedScore[bestPreviousSegmentIndex]:
                bestPreviousSegmentIndex = i
        result = [segments[bestPreviousSegmentIndex]]
        while (bestPreviousSegmentIndex := previousSegmentIndexes[bestPreviousSegmentIndex]) is not None:
            result.insert(0, segments[bestPreviousSegmentIndex])
        return result

    @staticmethod
    def __consecutivenessScore(previousSegment: AlignmentSegment, currentSegment: AlignmentSegment):
        queryLength = min(currentSegment.endPosition.query.siteId - currentSegment.startPosition.query.siteId,
                          abs(previousSegment.endPosition.query.siteId - previousSegment.startPosition.query.siteId))
        referenceDistance = currentSegment.startPosition.reference.siteId - previousSegment.endPosition.reference.siteId
        referenceLength = min(
            currentSegment.endPosition.reference.siteId - currentSegment.startPosition.reference.siteId,
            previousSegment.endPosition.reference.siteId - previousSegment.startPosition.reference.siteId)
        queryDistance = abs(currentSegment.startPosition.query.siteId - previousSegment.endPosition.query.siteId)
        if min(referenceLength + 2 * referenceDistance, queryLength + 2 * queryDistance) < 0:
            return -math.inf

        distanceSum = referenceDistance + queryDistance
        distanceDiff = referenceDistance - queryDistance
        return -(distanceSum ** 2 + distanceDiff ** 2) / max(abs(distanceSum), abs(distanceDiff), 1)
