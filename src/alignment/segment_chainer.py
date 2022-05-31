from __future__ import annotations

import math
from typing import Iterable, List

from src.alignment.segments import AlignmentSegment


class SegmentChainer:
    def chain(self, segments: Iterable[AlignmentSegment]):
        def initialOrderingKey(segment: AlignmentSegment):
            return segment.startPosition.reference.siteId + segment.endPosition.reference.siteId \
                   + segment.startPosition.query.siteId + segment.endPosition.query.siteId

        segments = sorted(segments, key=initialOrderingKey)
        cumulatedScore = [-math.inf] * len(segments)
        listOfPreviousSegmentIndexes: List[int | None] = [None] * len(segments)
        bestPreviousSegmentIndex = 0
        for i, currentSegment in enumerate(segments):
            cumulatedScore[i] = 0
            for j, previousSegment in enumerate(segments[:i]):
                currentScore = cumulatedScore[j] + self.__consecutivenessScore(previousSegment,
                                                                               currentSegment)
                if currentScore > cumulatedScore[i]:
                    cumulatedScore[i] = currentScore
                    listOfPreviousSegmentIndexes[i] = j
            cumulatedScore[i] += currentSegment.segmentScore
            if cumulatedScore[i] > cumulatedScore[bestPreviousSegmentIndex]:
                bestPreviousSegmentIndex = i
        result = [segments[bestPreviousSegmentIndex]]
        while (bestPreviousSegmentIndex := listOfPreviousSegmentIndexes[bestPreviousSegmentIndex]) is not None:
            result.insert(0, segments[bestPreviousSegmentIndex])
        return result

    def __consecutivenessScore(self, prevSegment: AlignmentSegment, currSegment: AlignmentSegment):
        queryLength = min(currSegment.endPosition.query.siteId - currSegment.startPosition.query.siteId,
                          prevSegment.endPosition.query.siteId - prevSegment.startPosition.query.siteId)
        referenceDistance = currSegment.startPosition.reference.siteId - prevSegment.endPosition.reference.siteId
        referenceLength = min(currSegment.endPosition.reference.siteId - currSegment.startPosition.reference.siteId,
                              prevSegment.endPosition.reference.siteId - prevSegment.startPosition.reference.siteId)
        queryDistance = currSegment.startPosition.query.siteId - prevSegment.endPosition.query.siteId
        if min(referenceLength + 2 * referenceDistance, queryLength + 2 * queryDistance) < 0:
            return -math.inf

        distanceSum = referenceDistance + queryDistance
        distanceDiff = referenceDistance - queryDistance
        return -(distanceSum ** 2 + distanceDiff ** 2) / max(abs(distanceSum), abs(distanceDiff), 1)
