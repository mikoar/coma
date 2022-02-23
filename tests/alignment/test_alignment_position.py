from __future__ import annotations

from typing import Tuple, List

import pytest

from src.alignment.alignment_position import NotAlignedPosition, NotAlignedQueryPosition, NotAlignedReferencePosition, \
    AlignedPair
from src.correlation.optical_map import PositionWithSiteId


class TestAlignedPair(AlignedPair):
    def __init__(self, reference: int, query: int, queryShift: int = 0):
        super().__init__(PositionWithSiteId(reference, 0), PositionWithSiteId(query, 0), queryShift)

    @staticmethod
    def createAlignedPairs(pairs: List[Tuple[int, int] | Tuple[int, int, int]]):
        return list(map(lambda p: TestAlignedPair(p[0], p[1], p[2] if 2 < len(p) else 0), pairs))


@pytest.mark.parametrize("position", [
    (NotAlignedQueryPosition(PositionWithSiteId(123, 0), 0)),
    (NotAlignedReferencePosition(PositionWithSiteId(234, 0)))
])
def test_score_notAlignedPosition_returnsUnmatchedPenalty(position: NotAlignedPosition):
    unmatchedPenalty = -100
    assert position.getScoredPosition(123, 2, unmatchedPenalty).score == unmatchedPenalty


def test_score_alignedPair_perfectMatch():
    position = TestAlignedPair(1, 1, 0)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScoredPosition(perfectMatchScore, scoreMultiplier, -123).score == 200


def test_score_alignedPair_distancePenalty():
    position = TestAlignedPair(1, 1, 50)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScoredPosition(perfectMatchScore, scoreMultiplier, -123).score == 100


if __name__ == '__main__':
    pytest.main(args=[__file__])
