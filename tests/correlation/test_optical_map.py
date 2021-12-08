from mimetypes import init

import pytest

from src.alignment.aligner import Aligner
from src.correlation.optical_map import OpticalMap
from src.correlation.sequence_generator import SequenceGenerator


def test_refineAlignment_correctPeakPosition():
    reference = OpticalMap(1, 1000, [20, 100, 110, 300, 310, 330, 400, 1000])
    query = OpticalMap(2, 100, [0, 10, 30, 100])

    initialGenerator = SequenceGenerator(2, 2)
    refineGenerator = SequenceGenerator(1, 1)

    initialAlignment = query.getInitialAlignment(reference, initialGenerator)
    refinedAlignment = initialAlignment.refine(refineGenerator, 10)

    assert initialAlignment.maxPeak.positionInReference == 350
    assert refinedAlignment.maxPeak.positionInReference == 350


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

    assert initialCorrelation.maxPeak.positionInReference == 9977344
    assert refinedCorrelation.maxPeak.positionInReference == 9977071

    initialAligner = Aligner(2 * initialResolution * initialBlur)
    refineAligner = Aligner(2 * refineResolution * refineBlur)
    nonRefinedAlignmentResult = initialAligner.align(reference, query, initialCorrelation.maxPeak.positionInReference)
    refinedAlignmentResult = refineAligner.align(reference, query, refinedCorrelation.maxPeak.positionInReference)

    assert refinedAlignmentResult.alignedPairs == nonRefinedAlignmentResult.alignedPairs


if __name__ == '__main__':
    pytest.main(args=[__file__])
