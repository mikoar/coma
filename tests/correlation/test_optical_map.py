from typing import List

import numpy as np
import pytest

from src.correlation.optical_map import OpticalMap, toRelativeGenomicPositions, PositionWithSiteId
from src.correlation.sequence_generator import SequenceGenerator


def test_trim():
    opticalMap = OpticalMap(2, 130, [20, 50, 100])
    trimmed = opticalMap.trim()
    assert trimmed.length == 81
    assert trimmed.positions == [0, 30, 80]


def test_trim_empty():
    opticalMap = OpticalMap(2, 20, [])
    trimmed = opticalMap.trim()
    assert trimmed.positions == []


def test_getPositionsWithSiteIds():
    opticalMap = OpticalMap(2, 101, [0, 10, 30, 100])
    assert list(opticalMap.getPositionsWithSiteIds()) \
           == [PositionWithSiteId(i, p) for i, p in [(1, 0), (2, 10), (3, 30), (4, 100)]]


def test_getPositionsWithSiteIds_reverse():
    opticalMap = OpticalMap(2, 101, [0, 10, 30, 100])
    assert list(opticalMap.getPositionsWithSiteIds(True)) \
           == [PositionWithSiteId(i, p) for i, p in [(4, 0), (3, 70), (2, 90), (1, 100)]]


@pytest.mark.parametrize("reverse", [True, False])
def test_getPositionsWithSiteIds_empty(reverse: bool):
    opticalMap = OpticalMap(2, 120, [])
    assert list(opticalMap.getPositionsWithSiteIds(reverse)) == []


def test_refineAlignment_correctPeakPosition():
    reference = OpticalMap(1, 1000, [20, 100, 110, 300, 310, 330, 400, 1000])
    query = OpticalMap(2, 100, [0, 10, 30, 100])

    initialGenerator = SequenceGenerator(2, 2)
    refineGenerator = SequenceGenerator(1, 1)

    initialAlignment = query.getInitialAlignment(reference, initialGenerator)
    refinedAlignment = initialAlignment.refine(refineGenerator, 10, 1)

    assert initialAlignment.maxPeak.position == 300
    assert refinedAlignment.maxPeak.position == 300


@pytest.mark.parametrize("positions,resolution,expected,start", [
    ([1, 3, 5], 1, [1, 3, 5], 0),
    ([1, 3, 5], 2, [2, 6, 10], 0),
    ([1, 3, 5], 3, [4, 10, 16], 0),
    ([1, 3, 5], 9, [13, 31, 49], 0),
    ([1, 3, 5], 10, [14, 34, 54], 0),
    ([1, 3, 5], 10, [34, 54, 74], 20),
])
def test_toRelativeGenomicPositions(positions: List[int], resolution: int, expected: List[int], start: int):
    assert toRelativeGenomicPositions(np.array(positions), resolution, start).tolist() == expected


def test_getInitialAlignment_whenQueryIsShorterThanReference_returnsEmptyResult():
    query = OpticalMap(1, 1000, [20, 100, 110, 300, 310, 330, 400, 1000])
    reference = OpticalMap(2, 100, [0, 10, 30, 100])

    initialAlignment = query.getInitialAlignment(reference, SequenceGenerator(2, 2))

    assert len(list(initialAlignment.peaks)) == 0
    assert initialAlignment.getScore() == 0


if __name__ == '__main__':
    pytest.main(args=[__file__])
