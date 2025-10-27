import itertools
from typing import List, Iterable

from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import PositionWithSiteId


class AlignmentSegmentConflictResolver:
    def __init__(self, segmentChainer: SegmentChainer):
        self.segmentChainer = segmentChainer

    def resolveConflicts(self, segments: List[AlignmentSegment], queryPositions: list[PositionWithSiteId]):
        if len(segments) < 2:
            return AlignmentSegmentsWithResolvedConflicts(segments)

        resolvedSegments = self.__pairAndResolveConflicts(segments, queryPositions)
        return AlignmentSegmentsWithResolvedConflicts(resolvedSegments)

    def __pairAndResolveConflicts(self, segments: Iterable[AlignmentSegment], queryPositions: list[PositionWithSiteId]):
        chainedSegments = self.segmentChainer.chain(segments, queryPositions)
        for (i0, i1) in self.__pairIndexes(len(chainedSegments)):
            pair = chainedSegments[i0].checkForConflicts(chainedSegments[i1])
            chainedSegments[i0], chainedSegments[i1] = pair.resolveConflict()
        return chainedSegments

    @staticmethod
    def __pairIndexes(length: int):
        a, b = itertools.tee(range(length))
        next(b, None)
        return zip(a, b)


class AlignmentSegmentsWithResolvedConflicts:
    def __init__(self, segments: List[AlignmentSegment]):
        self.segments = segments
