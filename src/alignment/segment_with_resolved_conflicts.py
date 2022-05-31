import itertools
from typing import List, Iterable

from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segments import AlignmentSegment


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
        segments = SegmentChainer().chain(segments)
        for (i0, i1) in AlignmentSegmentsWithResolvedConflicts.__pairIndexes(len(segments)):
            pair = segments[i0].checkForConflicts(segments[i1])
            segments[i0], segments[i1] = pair.resolveConflict()
        return segments

    @staticmethod
    def __pairIndexes(length: int):
        a, b = itertools.tee(range(length))
        next(b, None)
        return zip(a, b)
