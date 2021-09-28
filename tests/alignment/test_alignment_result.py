from __future__ import annotations

import pytest
from mock import Mock, call

from src.alignment.aligned_pair import AlignedPair
from src.alignment.alignment_result import AlignmentResult
from src.alignment.region_score_penalties import RegionScorePenalties


def test_getRegionScores_appliesPenalties_1pair():
    pairs = [AlignedPair(1, 1)]
    alignmentResult = AlignmentResult(1, 1, pairs)
    penalties = __getPenaltiesMock()
    penalties.getDistancePenalty = Mock(return_value=2)
    penalties.getUnmatchedLabelPenalty = Mock(return_value=5)

    perfectMatchScore = 100
    regionScores = alignmentResult.getRegionScores(penalties, perfectMatchScore)

    assert len(regionScores.scores) == 1
    assert regionScores.scores == [100 - 2 - 5]
    penalties.getDistancePenalty.assert_called_once_with(None, pairs[0])
    penalties.getUnmatchedLabelPenalty.assert_called_once_with(None, pairs[0])


def test_getRegionScores_appliesPenalties_2pairs():
    pairs = [AlignedPair(1, 1), AlignedPair(1, 1)]
    alignmentResult = AlignmentResult(1, 1, pairs)
    penalties = __getPenaltiesMock()
    penalties.getDistancePenalty = Mock(return_value=2)
    penalties.getUnmatchedLabelPenalty = Mock(return_value=5)

    perfectMatchScore = 100
    regionScores = alignmentResult.getRegionScores(penalties, perfectMatchScore)

    assert len(regionScores.scores) == 2
    assert regionScores.scores == [100 - 2 - 5] * 2
    calls = [call(None, pairs[0]), call(pairs[0], pairs[1])]
    penalties.getDistancePenalty.assert_has_calls(calls)
    penalties.getUnmatchedLabelPenalty.assert_has_calls(calls)


def __getPenaltiesMock() -> RegionScorePenalties | Mock:
    return Mock(spec=RegionScorePenalties)


if __name__ == '__main__':
    pytest.main(args=[__file__])
