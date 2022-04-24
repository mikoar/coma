from __future__ import annotations

import itertools
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
    def referenceStart(self):
        return self.alignedPositions[0].reference.siteId

    @property
    def referenceEnd(self):
        return self.alignedPositions[-1].reference.siteId + 1

    def checkForConflicts(self, other: AlignmentSegment):
        if self.endOverlapsWithStartOf(other):
            return _SegmentPairWithReferenceConflict.create(self, other, other.referenceStart,
                                                            self.referenceEnd)
        if other.endOverlapsWithStartOf(self):
            return _SegmentPairWithReferenceConflict.create(self, other, self.referenceStart,
                                                            other.referenceEnd)
        return _SegmentPairWithNoConflict(self, other)

    def endOverlapsWithStartOf(self, other):
        return other.referenceEnd >= self.referenceEnd > other.referenceStart >= self.referenceStart

    def sliceByReference(self, startSiteId: int, endSiteId) -> AlignmentSegment:
        slicedAtStart = itertools.dropwhile(
            lambda p: not isinstance(p, AlignedPair) or p.reference.siteId < startSiteId,
            self.positions)
        positions = list(itertools.takewhile(lambda p: not isinstance(p, AlignedPair) or p.reference.siteId < endSiteId,
                                             slicedAtStart))
        self.__trimNotAlignedPositionsFromEnd(positions)
        return AlignmentSegment(positions, sum(p.score for p in positions))

    @staticmethod
    def __trimNotAlignedPositionsFromEnd(positions):
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
    def referenceStart(self):
        return 0

    @property
    def referenceEnd(self):
        return 0

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


class _SegmentPairWithReferenceConflict(_SegmentPair):
    def __init__(self, segment1: AlignmentSegment, conflictingSubsegment1: AlignmentSegment,
                 segment2: AlignmentSegment, conflictingSubsegment2: AlignmentSegment):
        super().__init__(segment1, segment2)
        self.conflictingSubsegment1 = conflictingSubsegment1
        self.conflictingSubsegment2 = conflictingSubsegment2

    @staticmethod
    def create(segment1: AlignmentSegment, segment2: AlignmentSegment,
               conflictStartSiteId: int, conflictEndSiteId: int):
        conflictingSubsegment1 = segment1.sliceByReference(conflictStartSiteId, conflictEndSiteId)
        conflictingSubsegment2 = segment2.sliceByReference(conflictStartSiteId, conflictEndSiteId)
        return _SegmentPairWithReferenceConflict(segment1, conflictingSubsegment1, segment2, conflictingSubsegment2)

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if self.conflictingSubsegment1.segmentScore > self.conflictingSubsegment2.segmentScore:
            return self.segment1, self.segment2 - self.conflictingSubsegment2
        else:
            return self.segment1 - self.conflictingSubsegment1, self.segment2


class AlignmentSegmentsWithResolvedConflicts:
    @staticmethod
    def create(segments: List[AlignmentSegment]):
        if len(segments) < 2:
            return AlignmentSegmentsWithResolvedConflicts(segments)

        resolvedSegments = AlignmentSegmentsWithResolvedConflicts.__pairAndResolveConflicts(segments)
        notEmptySegments = [s for s in resolvedSegments if s != AlignmentSegment.empty]
        return AlignmentSegmentsWithResolvedConflicts(notEmptySegments)

    def __init__(self, segments: List[AlignmentSegment]):
        self.segments = segments

    @staticmethod
    def __pairAndResolveConflicts(segments: Iterable[AlignmentSegment]):
        segments = list(segments)
        segmentIndexPairs = itertools.combinations(range(len(segments)), 2)
        for (i0, i1) in segmentIndexPairs:
            pair = segments[i0].checkForConflicts(segments[i1])
            segments[i0], segments[i1] = pair.resolveConflict()
        return segments
