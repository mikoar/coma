from __future__ import annotations

import itertools
import math
from abc import ABC, abstractmethod
from typing import List, Tuple, Iterable

from src.alignment.alignment_position import ScoredAlignmentPosition, ScoredAlignedPair, AlignedPair


class AlignmentSegment:
    empty: AlignmentSegment

    def __init__(self, positions: List[ScoredAlignmentPosition], segmentScore: float):
        self.positions = positions
        self.segmentScore = segmentScore  # todo: refactor this
        self.alignedPositions = [p for p in positions if isinstance(p, ScoredAlignedPair)]

    @property
    def startPosition(self):
        return self.alignedPositions[0]

    @property
    def endPosition(self):
        return self.alignedPositions[-1]

    def checkForConflicts(self, other: AlignmentSegment):
        if self.endOverlapsWithStartOf(other):
            return _SegmentPairWithConflict.create(self, other, other.startPosition, self.endPosition)
        return _SegmentPairWithNoConflict(self, other)

    def endOverlapsWithStartOf(self, other: AlignmentSegment):
        return self.startPosition.lessOrEqualOnAnySequence(other.startPosition) or \
               other.startPosition.lessOrEqualOnAnySequence(self.endPosition) or \
               self.endPosition.lessOrEqualOnAnySequence(other.endPosition)

    def slice(self, start: AlignedPair, end: AlignedPair) -> AlignmentSegment:
        slicedAtStart = itertools.dropwhile(
            lambda p: not isinstance(p, AlignedPair) or p.lessOnBothSequences(start), self.positions)
        positions = list(
            itertools.takewhile(
                lambda p: not isinstance(p, AlignedPair) or p.lessOrEqualOnAnySequence(end),
                slicedAtStart))
        self.__trimNotAlignedPositionsFromEnd(positions)
        return AlignmentSegment(positions, sum(p.score for p in positions))

    @staticmethod
    def chain(segments: Iterable[AlignmentSegment]):
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
                currentScore = cumulatedScore[j] + AlignmentSegment.__consecutivenessScore(previousSegment,
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

    @staticmethod
    def __consecutivenessScore(prevSegment: AlignmentSegment, currSegment: AlignmentSegment):
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

    @staticmethod
    def __trimNotAlignedPositionsFromEnd(positions):
        if positions:
            while not isinstance(positions[-1], AlignedPair):
                positions.pop()

    def __eq__(self, other):
        return isinstance(other, AlignmentSegment) \
               and other.segmentScore == self.segmentScore \
               and other.positions == self.positions

    def __sub__(self, other: AlignmentSegment):
        positions = [p for p in self.positions if p not in other.positions]
        if not positions:
            return AlignmentSegment.empty
        return AlignmentSegment(positions, sum(p.score for p in positions))

    def __repr__(self):
        return f"score: {self.segmentScore}, positions: {self.positions}"


class __EmptyAlignmentSegment(AlignmentSegment):
    def __init__(self):
        super().__init__([], 0.)

    @property
    def startPosition(self):
        return AlignedPair.null

    @property
    def endPosition(self):
        return AlignedPair.null

    def checkForConflicts(self, other: AlignmentSegment):
        return _SegmentPairWithNoConflict(self, other)


AlignmentSegment.empty = __EmptyAlignmentSegment()


class _SegmentPair(ABC):
    def __init__(self, segment1: AlignmentSegment, segment2: AlignmentSegment):
        self.segment1 = segment1
        self.segment2 = segment2

    @abstractmethod
    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        pass


class _SegmentPairWithNoConflict(_SegmentPair):
    def __init__(self, segment1: AlignmentSegment, segment2: AlignmentSegment):
        super().__init__(segment1, segment2)

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        return self.segment1, self.segment2


class _SegmentPairWithConflict(_SegmentPair):
    def __init__(self, segment1: AlignmentSegment, conflictingSubsegment1: AlignmentSegment,
                 segment2: AlignmentSegment, conflictingSubsegment2: AlignmentSegment):
        super().__init__(segment1, segment2)
        self.conflictingSubsegment1 = conflictingSubsegment1
        self.conflictingSubsegment2 = conflictingSubsegment2

    @staticmethod
    def create(segment1: AlignmentSegment, segment2: AlignmentSegment,
               conflictStart: AlignedPair, conflictEnd: AlignedPair):
        conflictingSubsegment1 = segment1.slice(conflictStart, conflictEnd)
        conflictingSubsegment2 = segment2.slice(conflictStart, conflictEnd)
        return _SegmentPairWithConflict(segment1, conflictingSubsegment1, segment2, conflictingSubsegment2)

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if self.conflictingSubsegment1.segmentScore > self.conflictingSubsegment2.segmentScore:
            return self.segment1, self.segment2 - self.conflictingSubsegment2
        else:
            return self.segment1 - self.conflictingSubsegment1, self.segment2
