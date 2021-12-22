from typing import List

import numpy as np
import pytest

from src.alignment.aligner import Aligner
from src.correlation.optical_map import OpticalMap, adjustPeakPositions
from src.correlation.sequence_generator import SequenceGenerator


def test_refineAlignment_correctPeakPosition():
    reference = OpticalMap(1, 1000, [20, 100, 110, 300, 310, 330, 400, 1000])
    query = OpticalMap(2, 100, [0, 10, 30, 100])

    initialGenerator = SequenceGenerator(2, 2)
    refineGenerator = SequenceGenerator(1, 1)

    initialAlignment = query.getInitialAlignment(reference, initialGenerator)
    refinedAlignment = initialAlignment.refine(refineGenerator, 10)

    assert initialAlignment.maxPeak.position == 350
    assert refinedAlignment.maxPeak.position == 350


@pytest.mark.skip()
def test_refineAlignment_snapshot_realData():
    reference = OpticalMap(21, 10102331,
                           [9922485, 9926013, 9942837, 9958552, 9962567, 9980530, 9981296, 10008776, 10017980,
                            10026400,
                            10030478, 10102331])
    query = OpticalMap(1955058, 107027, [1862, 18773, 34673, 38770, 56914, 77567, 84923, 94390, 102845, 107027])
    initialResolution = 256
    initialBlur = 4
    initialGenerator = SequenceGenerator(initialResolution, initialBlur)
    refineResolution = 40
    refineBlur = 2
    refineGenerator = SequenceGenerator(refineResolution, refineBlur)
    initialCorrelation = query.getInitialAlignment(reference, initialGenerator)
    refinedCorrelation = initialCorrelation.refine(refineGenerator, 3000)

    assert initialCorrelation.maxPeak.position == 9977344
    assert refinedCorrelation.maxPeak.position == 9977071

    initialAligner = Aligner(2 * initialResolution * initialBlur)
    refineAligner = Aligner(2 * refineResolution * refineBlur)
    nonRefinedAlignmentResult = initialAligner.align(reference, query, initialCorrelation.maxPeak.position)
    refinedAlignmentResult = refineAligner.align(reference, query, refinedCorrelation.maxPeak.position)

    assert refinedAlignmentResult.alignedPairs == nonRefinedAlignmentResult.alignedPairs


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


if __name__ == '__main__':
    pytest.main(args=[__file__])
