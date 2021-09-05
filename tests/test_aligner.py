
import pytest
from src.aligner import Aligner
from src.optical_map import OpticalMap


@pytest.mark.parametrize("reference,query,peakPosition", [
    (
        OpticalMap(1, 9, [0, 4, 8]),
        OpticalMap(1, 9, [0, 4, 8]),
        4
    ),
    (
        OpticalMap(1, 10, [0, 5, 9]),
        OpticalMap(1, 10, [0, 5, 9]),
        5
    ),
    (
        OpticalMap(1, 30, [10, 15, 19]),
        OpticalMap(1, 10, [0, 5, 9]),
        15
    )
])
def test_simple(reference, query, peakPosition):
    result = Aligner(0).align(reference, query, peakPosition)

    assert [(1, 1), (2, 2), (3, 3)] == result.alignedPairs


def test_simple_reverseStrand():
    reference = OpticalMap(1, 9, [0, 4, 8])
    query = OpticalMap(1, 9, [0, 4, 8])

    result = Aligner(0).align(reference, query, 4, True)

    assert [(1, 3), (2, 2), (3, 1)] == result.alignedPairs


def test_extraPositionsOnReference():
    reference = OpticalMap(1, 30, [9, 10, 15, 19, 20])
    query = OpticalMap(1, 10, [0, 5, 9])

    result = Aligner(0).align(reference, query, 15)

    assert 10 == result.referenceStartPosition
    assert [(2, 1), (3, 2), (4, 3)] == result.alignedPairs


def test_alignWithinMaxDistance():
    reference = OpticalMap(1, 300, [88, 89, 149, 198, 199])
    query = OpticalMap(1, 100, [0, 49, 89])

    result = Aligner(10).align(reference, query, 149)

    assert [(2, 1), (3, 2), (4, 3)] == result.alignedPairs


def test_lastPositionAboveMaxDistance():
    reference = OpticalMap(1, 300, [99, 149, 199])
    query = OpticalMap(1, 100, [0, 49, 89])

    result = Aligner(10).align(reference, query, 149)

    assert 199 == result.referenceEndPosition
    assert [(1, 1), (2, 2)] == result.alignedPairs


def test_deletion():
    reference = OpticalMap(1, 300, [99, 149, 174, 189])
    query = OpticalMap(1, 100, [0, 49, 89])

    result = Aligner(10).align(reference, query, 149)

    assert 199 == result.referenceEndPosition
    assert [(1, 1), (2, 2), (4, 3)] == result.alignedPairs


if __name__ == '__main__':
    pytest.main(args=[__file__])
