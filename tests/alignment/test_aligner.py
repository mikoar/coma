from typing import Callable

import pytest

from src.alignment.aligned_pair import AlignedPair
from src.alignment.aligner import Aligner
from src.correlation.optical_map import OpticalMap


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

    assert result.alignedPairs == [(1, 1), (2, 2), (3, 3)]


def test_perfectMatch_reverseStrand():
    reference = OpticalMap(1, length=9, positions=[0, 4, 8])
    query = OpticalMap(1, length=9, positions=[0, 4, 8])

    result = Aligner(0).align(reference, query, 4, True)

    assert result.alignedPairs == [(1, 3), (2, 2), (3, 1)]


def test_ignoresExtraPositionsOnReferenceBeforeAndAfterAlignment():
    reference = OpticalMap(1, length=30, positions=[9, 10, 15, 19, 20])
    query = OpticalMap(1, length=10, positions=[0, 5, 9])

    result = Aligner(0).align(reference, query, 15)

    assert result.referenceStartPosition == 10
    assert result.alignedPairs == [(2, 1), (3, 2), (4, 3)]


def test_alignsPositionsWithinMaxDistanceOnly():
    maxDistance = 10
    reference = OpticalMap(1, length=300, positions=[89, 115, 150, 185, 210])
    query = OpticalMap(1, length=100, positions=[0, 25, 50, 75, 99])

    result = Aligner(maxDistance).align(reference, query, 150)

    assert result.alignedPairs == [(2, 2), (3, 3), (4, 4)]


def test_returnsCorrectQueryShifts():
    maxDistance = 10
    reference = OpticalMap(1, length=300, positions=[105, 115, 150, 185, 194])
    query = OpticalMap(1, length=100, positions=[0, 25, 50, 75, 99])

    result = Aligner(maxDistance).align(reference, query, 150)

    distanceSelector: Callable[[AlignedPair], int] = lambda pair: pair.queryShift
    assert list(map(distanceSelector, result.alignedPairs)) == [-5, 10, 0, -10, 5]


def test_ignoresPositionBeyondMaxDistance():
    reference = OpticalMap(1, length=300, positions=[99, 149, 199])
    query = OpticalMap(1, length=100, positions=[0, 49, 89])
    maxDistance = 10

    result = Aligner(maxDistance).align(reference, query, 149)

    assert result.referenceEndPosition == 199
    assert result.alignedPairs == [(1, 1), (2, 2)]


@pytest.mark.parametrize("reference,query", [
    (
            OpticalMap(1, length=300, positions=[100, 149, 174, 189]),
            OpticalMap(1, length=100, positions=[0, 49, 89])
    )
])
def test_handlesDeletions_referencePosition3OutOfRange_NotAligned(reference, query):
    result = Aligner(10).align(reference, query, 150)

    assert result.referenceEndPosition == 200
    assert result.alignedPairs == [(1, 1), (2, 2), (4, 3)]


@pytest.mark.parametrize("reference,query,expected", [
    (
            OpticalMap(1, length=300, positions=[100, 149, 150, 189]),
            OpticalMap(1, length=100, positions=[0, 49, 89]),
            [(1, 1), (2, 2), (3, 2), (4, 3)]
    ), (
            OpticalMap(1, length=300, positions=[100, 149, 188, 189]),
            OpticalMap(1, length=100, positions=[0, 49, 89]),
            [(1, 1), (2, 2), (3, 3), (4, 3)]
    )
])
def test_handlesDeletions_referencePosition3InRange_Aligns2ReferencesTo1Query(reference, query, expected):
    result = Aligner(10).align(reference, query, 150)

    assert result.referenceEndPosition == 200
    assert result.alignedPairs == expected


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

    assert result.referenceEndPosition == 200
    assert result.alignedPairs == [(1, 1), (2, 2), (3, 4)]


if __name__ == '__main__':
    pytest.main(args=[__file__])
