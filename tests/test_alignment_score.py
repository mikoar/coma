from __future__ import annotations
from typing import List, Union
from mock.mock import call
import pytest
from mock import Mock
from src.alignment_score import AlignmentScorer, RegionScorePenalties, RegionScorer, AlignmentResult, AlignedPair


def test_fullAlignment():
    alignmentScores = [5., 6., 7., 4., 3.]
    segment = AlignmentScorer().getSegmentWithMaxScore(alignmentScores)
    assert 25. == segment.score
    assert 0 == segment.start
    assert 4 == segment.end


def test_partialAlignment():
    alignmentScores = [-1., 5., 6., 7., 4., 3., -1]
    segment = AlignmentScorer().getSegmentWithMaxScore(alignmentScores)
    assert 25. == segment.score
    assert 1 == segment.start
    assert 5 == segment.end


def test_alignmentWithGap():
    alignmentScores = [1., 1., -3., 2., 1., -3., 2.]
    segment = AlignmentScorer().getSegmentWithMaxScore(alignmentScores)
    assert 3. == segment.score
    assert 3 == segment.start
    assert 4 == segment.end


def test_getRegionScores_appliesPenalties_1pair():
    pairs = [AlignedPair(1, 1)]
    alignmentResult = AlignmentResult(1, 1, pairs)
    penalties = __getPenaltiesMock()
    penalties.getDistancePenalty = Mock(return_value=2)
    penalties.getUnmatchedLabelPenalty = Mock(return_value=5)

    perfectMatchScore = 100
    regionScores = list(RegionScorer(penalties, perfectMatchScore).getRegionScores(alignmentResult))

    assert len(regionScores) == 1
    assert regionScores == [100 - 2 - 5]
    penalties.getDistancePenalty.assert_called_once_with(None, pairs[0])
    penalties.getUnmatchedLabelPenalty.assert_called_once_with(None, pairs[0])


def test_getRegionScores_appliesPenalties_2pairs():
    pairs = [AlignedPair(1, 1), AlignedPair(1, 1)]
    alignmentResult = AlignmentResult(1, 1, pairs)
    penalties = __getPenaltiesMock()
    penalties.getDistancePenalty = Mock(return_value=2)
    penalties.getUnmatchedLabelPenalty = Mock(return_value=5)

    perfectMatchScore = 100
    regionScores = list(RegionScorer(penalties, perfectMatchScore).getRegionScores(alignmentResult))

    assert len(regionScores) == 2
    assert regionScores == [100 - 2 - 5] * 2
    calls = [call(None, pairs[0]), call(pairs[0], pairs[1])]
    penalties.getDistancePenalty.assert_has_calls(calls)
    penalties.getUnmatchedLabelPenalty.assert_has_calls(calls)


@pytest.mark.parametrize("previousPair, pair, unmatchedCount", [
    (None, AlignedPair(1, 1), 0),
    (None, AlignedPair(30, 1), 0),
    (AlignedPair(1, 1), AlignedPair(2, 2), 0),
    (AlignedPair(100, 10), AlignedPair(101, 11), 0),
    (AlignedPair(1, 1), AlignedPair(2, 3), 1),
    (AlignedPair(1, 1), AlignedPair(3, 2), 1),
    (AlignedPair(100, 10), AlignedPair(101, 12), 1),
    (AlignedPair(1, 1), AlignedPair(3, 3), 2),
    (AlignedPair(1, 3), AlignedPair(2, 2), 0),
])
def test_regionScorePenalties_unmatchedLabelPenalty(previousPair: AlignedPair | None, pair: AlignedPair, unmatchedCount: int):
    unmatchedLabelPenalty = 100
    penalties = RegionScorePenalties(unmatchedLabelPenalty)
    assert penalties.getUnmatchedLabelPenalty(previousPair, pair) == unmatchedCount * unmatchedLabelPenalty


# def test_regionScorePenalties_distancePenalty():
#     alignmentResult = AlignmentResult(1, 1000, [
#         AlignedPair(1, 1),
#         AlignedPair(2, 2),
#         AlignedPair(3, 3)
#     ])
#     perfectMatchScore = 123
#     distancePenaltyExponent = 2
#     regionScores = RegionScorer(perfectMatchScore, distancePenaltyExponent).getRegionScores(alignmentResult)

#     assert list(regionScores) == [perfectMatchScore, perfectMatchScore, perfectMatchScore]


def __getPenaltiesMock() -> Union[RegionScorePenalties, Mock]:
    return Mock(spec=RegionScorePenalties)


if __name__ == '__main__':
    pytest.main(args=[__file__])
