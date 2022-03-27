from __future__ import annotations

import pytest

from src.alignment.alignment_position import NotAlignedPosition, NotAlignedQueryPosition, NotAlignedReferencePosition
from src.correlation.optical_map import PositionWithSiteId
from tests.test_doubles.alignment_segment_stub import AlignedPairStub


@pytest.mark.parametrize("position", [
    (NotAlignedQueryPosition(PositionWithSiteId(123, 0), 0)),
    (NotAlignedReferencePosition(PositionWithSiteId(234, 0)))
])
def test_score_notAlignedPosition_returnsUnmatchedPenalty(position: NotAlignedPosition):
    unmatchedPenalty = -100
    assert position.getScoredPosition(123, 2, unmatchedPenalty).score == unmatchedPenalty


def test_score_alignedPair_perfectMatch():
    position = AlignedPairStub(1, 1, 0)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScoredPosition(perfectMatchScore, scoreMultiplier, -123).score == 200


def test_score_alignedPair_distancePenalty():
    position = AlignedPairStub(1, 1, 50)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScoredPosition(perfectMatchScore, scoreMultiplier, -123).score == 100


if __name__ == '__main__':
    pytest.main(args=[__file__])
