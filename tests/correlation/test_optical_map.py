from typing import List

import numpy as np
import pytest

from src.correlation.optical_map import OpticalMap, adjustPeakPositions
from src.correlation.sequence_generator import SequenceGenerator


def test_refineAlignment_correctPeakPosition():
    reference = OpticalMap(1, 1000, [20, 100, 110, 300, 310, 330, 400, 1000])
    query = OpticalMap(2, 100, [0, 10, 30, 100])

    initialGenerator = SequenceGenerator(2, 2)
    refineGenerator = SequenceGenerator(1, 1)

    initialAlignment = query.getInitialAlignment(reference, initialGenerator)
    refinedAlignment = initialAlignment.refine(refineGenerator, 10, 1)

    assert initialAlignment.maxPeak.position == 300
    assert refinedAlignment.maxPeak.position == 300


@pytest.mark.parametrize("positions,resolution,adjustedPositions,start", [
    ([1, 3, 5], 1, [1, 3, 5], 0),
    ([1, 3, 5], 2, [2, 6, 10], 0),
    ([1, 3, 5], 3, [4, 10, 16], 0),
    ([1, 3, 5], 9, [13, 31, 49], 0),
    ([1, 3, 5], 10, [14, 34, 54], 0),
    ([1, 3, 5], 10, [34, 54, 74], 20),
])
def test_adjustPeakPositions(positions: List[int], resolution: int, adjustedPositions: List[int], start: int):
    assert adjustPeakPositions(np.array(positions), resolution, start).tolist() == adjustedPositions


def test_getInitialAlignment_whenQueryIsShorterThanReference_returnsEmptyResult():
    query = OpticalMap(1, 1000, [20, 100, 110, 300, 310, 330, 400, 1000])
    reference = OpticalMap(2, 100, [0, 10, 30, 100])

    initialAlignment = query.getInitialAlignment(reference, SequenceGenerator(2, 2))

    assert len(list(initialAlignment.peaks)) == 0
    assert initialAlignment.getScore() == 0


if __name__ == '__main__':
    pytest.main(args=[__file__])
