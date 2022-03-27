from __future__ import annotations

from typing import List, Tuple

from src.alignment.alignment_position import ScoredAlignedPair, AlignedPair
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import PositionWithSiteId


class AlignmentSegmentStub(AlignmentSegment):
    def __init__(self, positions: List[ScoredAlignedPairStub]):
        super().__init__(positions, 0, 0, 0)

    @staticmethod
    def createFromPairs(pairs: List[Tuple[int, int] | Tuple[int, int, int]]):
        return [AlignmentSegmentStub(
            list(map(lambda p: ScoredAlignedPairStub(p[0], p[1], p[2] if 2 < len(p) else 0), pairs)))]


class ScoredAlignedPairStub(ScoredAlignedPair):
    def __init__(self, reference: int, query: int, queryShift: int = 0):
        super().__init__(AlignedPairStub(reference, query, queryShift), 0.)


class AlignedPairStub(AlignedPair):
    def __init__(self, reference: int, query: int, queryShift: int = 0):
        super().__init__(PositionWithSiteId(reference, 0), PositionWithSiteId(query, 0), queryShift)
