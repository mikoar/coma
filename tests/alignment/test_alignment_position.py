import pytest

from src.alignment.alignment_position import NotAlignedPosition, NotAlignedQueryPosition, NotAlignedReferencePosition, \
    AlignedPair
from src.correlation.optical_map import PositionWithSiteId


class TestAlignedPair(AlignedPair):
    def __init__(self, reference: int, query: int, queryShift: int = 0):
        super().__init__(PositionWithSiteId(reference, 0), PositionWithSiteId(query, 0), queryShift)


@pytest.mark.parametrize("position", [
    (NotAlignedQueryPosition(PositionWithSiteId(123, 0), 0)), (NotAlignedReferencePosition(PositionWithSiteId(234, 0)))
])
def test_score_notAlignedPosition_returnsUnmatchedPenalty(position: NotAlignedPosition):
    unmatchedPenalty = -100
    assert position.getScore(123, 2, unmatchedPenalty) == unmatchedPenalty


def test_score_alignedPair_perfectMatch():
    position = TestAlignedPair(1, 1, 0)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScore(perfectMatchScore, scoreMultiplier, -123) == 200


def test_score_alignedPair_distancePenalty():
    position = TestAlignedPair(1, 1, 50)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScore(perfectMatchScore, scoreMultiplier, -123) == 100


if __name__ == '__main__':
    pytest.main(args=[__file__])
