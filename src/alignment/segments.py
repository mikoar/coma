from __future__ import annotations

import itertools
from abc import ABC, abstractmethod
from typing import List, Tuple

from src.alignment.alignment_position import ScoredAlignmentPosition, ScoredAlignedPair, AlignedPair
from src.correlation.peak import Peak


class AlignmentSegment:
    empty: AlignmentSegment

    def __init__(self, positions: List[ScoredAlignmentPosition], segmentScore: float, peak: Peak):
        self.positions = positions
        self.segmentScore = segmentScore
        self.alignedPositions = [p for p in positions if isinstance(p, ScoredAlignedPair)]
        self.peak = peak

    @property
    def startPosition(self):
        return self.alignedPositions[0]

    @property
    def endPosition(self):
        return self.alignedPositions[-1]

    @property
    def reverse(self):
        return self.startPosition.query.siteId > self.endPosition.query.siteId

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
        return AlignmentSegment(positions, sum(p.score for p in positions), self.peak)

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
        return AlignmentSegment(positions, sum(p.score for p in positions), self.peak)

    def __repr__(self):
        return f"score: {self.segmentScore}, positions: {self.positions}"


class __EmptyAlignmentSegment(AlignmentSegment):
    def __init__(self):
        super().__init__([], 0., Peak.null)

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
