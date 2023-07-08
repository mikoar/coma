from __future__ import annotations

from typing import List, Tuple

from src.alignment.alignment_position import ScoredAlignedPair, AlignedPair, ScoredNotAlignedPosition, \
    NotAlignedReferencePosition, NotAlignedQueryPosition
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import PositionWithSiteId
from src.correlation.peak import Peak


def createPositionWithSiteId(position: Tuple[int, int] | int):
    return PositionWithSiteId(position[0], position[1]) \
        if type(position) is tuple \
        else PositionWithSiteId(position, position * 100)


class AlignmentSegmentStub(AlignmentSegment):
    def __init__(self, positions: List[ScoredAlignedPairStub], score: float = 0.):
        super().__init__(positions, score, Peak.null, positions)

    @staticmethod
    def createFromPairs(scoredPositionTuples: List[Tuple[int | None, int | None]
                                                   | Tuple[Tuple[int, int] | int | None,
                                                           Tuple[int, int] | int | None, float]]):
        positions = list(map(
            lambda pair: ScoredNotAlignedPositionStub(pair[0], pair[1], pair[2] if len(pair) > 2 else 0)
            if not pair[0] or not pair[1] else
            ScoredAlignedPairStub(pair[0], pair[1], 0, pair[2] if len(pair) > 2 else 0),
            scoredPositionTuples))
        return AlignmentSegmentStub(positions, sum(p.score for p in positions))


class ScoredAlignedPairStub(ScoredAlignedPair):
    def __init__(self, reference: Tuple[int, int] | int, query: Tuple[int, int] | int, queryShift: int = 0,
                 score: float = 0.):
        super().__init__(AlignedPairStub(reference, query, queryShift), score)


class ScoredNotAlignedPositionStub(ScoredNotAlignedPosition):
    def __init__(self, reference: Tuple[int, int] | int = None, query: Tuple[int, int] | int = None, score: float = 0.):
        if reference:
            super().__init__(NotAlignedReferencePosition(createPositionWithSiteId(reference)), score)
        else:
            super().__init__(NotAlignedQueryPosition(createPositionWithSiteId(query), 0), score)


class AlignedPairStub(AlignedPair):
    def __init__(self, reference: Tuple[int, int] | int, query: Tuple[int, int] | int, queryShift: int = 0):
        super().__init__(createPositionWithSiteId(reference), createPositionWithSiteId(query),
                         queryShift)
