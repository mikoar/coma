from __future__ import annotations

from typing import List, Tuple

from src.alignment.alignment_position import ScoredAlignedPair, AlignedPair, ScoredNotAlignedPosition, \
    NotAlignedReferencePosition, NotAlignedQueryPosition
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import PositionWithSiteId


class AlignmentSegmentStub(AlignmentSegment):
    def __init__(self, positions: List[ScoredAlignedPairStub], score: float = 0.):
        super().__init__(positions, score)

    @staticmethod
    def createFromPairs(
            scoredPositionTuples: List[Tuple[int | None, int | None] | Tuple[int | None, int | None, float]]):
        positions = list(map(
            lambda pair: ScoredNotAlignedPositionStub(pair[0], pair[1], pair[2] if len(pair) > 2 else 0)
            if not pair[0] or not pair[1] else
            ScoredAlignedPairStub(pair[0], pair[1], 0, pair[2] if len(pair) > 2 else 0),
            scoredPositionTuples))
        return AlignmentSegmentStub(positions, sum(p.score for p in positions))


class ScoredAlignedPairStub(ScoredAlignedPair):
    def __init__(self, reference: int, query: int, queryShift: int = 0, score: float = 0.):
        super().__init__(AlignedPairStub(reference, query, queryShift), score)


class ScoredNotAlignedPositionStub(ScoredNotAlignedPosition):
    def __init__(self, reference: int = None, query: int = None, score: float = 0.):
        if reference:
            super().__init__(NotAlignedReferencePosition(PositionWithSiteId(reference, reference * 100)), score)
        else:
            super().__init__(NotAlignedQueryPosition(PositionWithSiteId(query, query * 100), 0), score)


class AlignedPairStub(AlignedPair):
    def __init__(self, reference: int, query: int, queryShift: int = 0):
        super().__init__(PositionWithSiteId(reference, reference * 100), PositionWithSiteId(query, query * 100),
                         queryShift)
