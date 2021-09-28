from __future__ import annotations

import pytest

from src.alignment.aligned_pair import AlignedPair
from src.alignment.region_score_penalties import RegionScorePenalties


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
def test_regionScorePenalties_unmatchedLabelPenalty(previousPair: AlignedPair | None, pair: AlignedPair,
                                                    unmatchedCount: int):
    unmatchedLabelPenalty = 100
    penalties = RegionScorePenalties(unmatchedLabelPenalty)
    assert penalties.getUnmatchedLabelPenalty(previousPair, pair) == unmatchedCount * unmatchedLabelPenalty


@pytest.mark.parametrize("previousPair, pair, distance", [
    (None, AlignedPair(1, 1, 0), 0),
    (None, AlignedPair(1, 1, 5), 5),
    (None, AlignedPair(1, 1, -5), 5),
    (AlignedPair(1, 1, 5), AlignedPair(1, 1, 5), 0),
    (AlignedPair(1, 1, -5), AlignedPair(1, 1, 5), 10),
])
def test_regionScorePenalties_distancePenalty(previousPair: AlignedPair | None, pair: AlignedPair,
                                              distance: int):
    penalties = RegionScorePenalties(distancePenaltyExponent=1)
    assert penalties.getDistancePenalty(previousPair, pair) == distance


@pytest.mark.parametrize("distancePenaltyExponent,penalty", [
    (2., 4.),
    (4., 16.),
])
def test_regionScorePenalties_distancePenalty_exponentialScaling(distancePenaltyExponent: float, penalty: float):
    penalties = RegionScorePenalties(distancePenaltyExponent=distancePenaltyExponent)
    assert penalties.getDistancePenalty(AlignedPair(1, 1, 0), AlignedPair(1, 1, 2)) == penalty


if __name__ == '__main__':
    pytest.main(args=[__file__])
