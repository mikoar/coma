import itertools
from typing import List, Iterable

from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segments import AlignmentSegment


class AlignmentSegmentConflictResolver:
    def __init__(self, segmentChainer: SegmentChainer):
        self.segmentChainer = segmentChainer

    def resolveConflicts(self, segments: List[AlignmentSegment]):
        if len(segments) < 2:
            return AlignmentSegmentsWithResolvedConflicts(segments)

        resolvedSegments = self.__pairAndResolveConflicts(segments)
        notEmptySegments = [s for s in resolvedSegments if s != AlignmentSegment.empty]
        return AlignmentSegmentsWithResolvedConflicts(notEmptySegments)

    def __pairAndResolveConflicts(self, segments: Iterable[AlignmentSegment]):
        chainedSegments = self.segmentChainer.chain(segments)
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
