
import pytest
from src.aligner import Aligner
from src.optical_map import OpticalMap


@pytest.mark.parametrize("reference,query,peakPosition", [
    (
        OpticalMap(1, length=9, positions=[0, 4, 8]),
        OpticalMap(1, length=9, positions=[0, 4, 8]),
        4
    ),
    (
        OpticalMap(1, length=10, positions=[0, 5, 9]),
        OpticalMap(1, length=10, positions=[0, 5, 9]),
        5
    ),
    (
        OpticalMap(1, length=30, positions=[10, 15, 19]),
        OpticalMap(1, length=10, positions=[0, 5, 9]),
        15
    )
])
def test_perfectMatch(reference, query, peakPosition):
    result = Aligner(0).align(reference, query, peakPosition)

    assert [(1, 1), (2, 2), (3, 3)] == result.alignedPairs


def test_perfectMatch_reverseStrand():
    reference = OpticalMap(1, length=9, positions=[0, 4, 8])
    query = OpticalMap(1, length=9, positions=[0, 4, 8])

    result = Aligner(0).align(reference, query, 4, True)

    assert [(1, 3), (2, 2), (3, 1)] == result.alignedPairs


def test_ignoresExtraPositionsOnReferenceBeforeAndAfterAlignment():
    reference = OpticalMap(1, length=30, positions=[9, 10, 15, 19, 20])
    query = OpticalMap(1, length=10, positions=[0, 5, 9])

    result = Aligner(0).align(reference, query, 15)

    assert 10 == result.referenceStartPosition
    assert [(2, 1), (3, 2), (4, 3)] == result.alignedPairs


def test_alignsPositionsWithinMaxDistanceOnly():
    reference = OpticalMap(1, length=300, positions=[88, 89, 149, 198, 199])
    query = OpticalMap(1, length=100, positions=[0, 49, 89])
    maxDistance = 10

    result = Aligner(maxDistance).align(reference, query, 149)

    assert [(2, 1), (3, 2), (4, 3)] == result.alignedPairs


def test_ignoresPositionBeyondMaxDistance():
    reference = OpticalMap(1, length=300, positions=[99, 149, 199])
    query = OpticalMap(1, length=100, positions=[0, 49, 89])
    maxDistance = 10

    result = Aligner(maxDistance).align(reference, query, 149)

    assert 199 == result.referenceEndPosition
    assert [(1, 1), (2, 2)] == result.alignedPairs


@pytest.mark.parametrize("reference,query", [
    (
        OpticalMap(1, length=300, positions=[100, 149, 174, 189]),
        OpticalMap(1, length=100, positions=[0, 49, 89])
    ), (
        OpticalMap(1, length=300, positions=[100, 149, 150, 189]),
        OpticalMap(1, length=100, positions=[0, 49, 89])
    ), (
        OpticalMap(1, length=300, positions=[100, 149, 188, 189]),
        OpticalMap(1, length=100, positions=[0, 49, 89])
    )
])
def test_handlesDeletions_noAlignmentForReferencePosition3(reference, query):
    result = Aligner(10).align(reference, query, 150)

    assert 200 == result.referenceEndPosition
    assert [(1, 1), (2, 2), (4, 3)] == result.alignedPairs


@pytest.mark.parametrize("reference,query", [
    (
        OpticalMap(1, length=300, positions=[100, 149, 189]),
        OpticalMap(1, length=100, positions=[0, 49, 74, 89])
    ), (
        OpticalMap(1, length=300, positions=[100, 149, 189]),
        OpticalMap(1, length=100, positions=[0, 49, 50, 89])
    ), (
        OpticalMap(1, length=300, positions=[100, 149, 189]),
        OpticalMap(1, length=100, positions=[0, 49, 88, 89])
    )
])
def test_handlesInsertions_noAlignmentForQueryPosition3(reference, query):
    result = Aligner(10).align(reference, query, 150)

    assert 200 == result.referenceEndPosition
    assert [(1, 1), (2, 2), (3, 4)] == result.alignedPairs


if __name__ == '__main__':
    pytest.main(args=[__file__])
