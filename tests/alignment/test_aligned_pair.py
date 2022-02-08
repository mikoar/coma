import pytest

from src.alignment.aligned_pair import NotAlignedPosition, NotAlignedQueryPosition, NotAlignedReferencePosition, \
    AlignedPair


@pytest.mark.parametrize("position", [
    (NotAlignedQueryPosition(123)), (NotAlignedReferencePosition(234))
])
def test_score_notAlignedPosition_returnsUnmatchedPenalty(position: NotAlignedPosition):
    unmatchedPenalty = -100
    assert position.getScore(123, 2, unmatchedPenalty) == unmatchedPenalty


def test_score_alignedPair_perfectMatch():
    position = AlignedPair(1, 1, 0)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScore(perfectMatchScore, scoreMultiplier, -123) == 200


def test_score_alignedPair_distancePenalty():
    position = AlignedPair(1, 1, 50)
    perfectMatchScore = 100
    scoreMultiplier = 2
    assert position.getScore(perfectMatchScore, scoreMultiplier, -123) == 100


if __name__ == '__main__':
    pytest.main(args=[__file__])
