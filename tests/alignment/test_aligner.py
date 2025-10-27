from typing import List
from unittest.mock import Mock

import pytest

from src.alignment.aligner import Aligner, AlignerEngine
from src.alignment.alignment_position import AlignedPair, AlignmentPosition
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver, \
    AlignmentSegmentsWithResolvedConflicts
from src.alignment.segments import AlignmentSegment
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.correlation.optical_map import OpticalMap
from src.correlation.peak import Peak


def getSut(maxDistance=0):
    segmentsFactoryMock: AlignmentSegmentsFactory = Mock(spec=AlignmentSegmentsFactory)
    segmentsFactoryMock.getSegments = lambda positions, peak: [AlignmentSegment(positions, 300, peak, positions)]
    segmentConflictResolverMock: AlignmentSegmentConflictResolver = Mock(spec=AlignmentSegmentConflictResolver)
    segmentConflictResolverMock.resolveConflicts = lambda segments, queryPosistions: AlignmentSegmentsWithResolvedConflicts(segments)
    return Aligner(AlignmentPositionScorer(100, 1, 0),
                   segmentsFactoryMock,
                   AlignerEngine(maxDistance),
                   segmentConflictResolverMock)


def __peak(position: int = 0):
    return Peak(position, 100)


@pytest.mark.parametrize("reference,query", [
    (
            OpticalMap(1, length=0, positions=[]),
            OpticalMap(1, length=0, positions=[])
    ),
    (
            OpticalMap(1, length=10, positions=[0, 5, 9]),
            OpticalMap(1, length=10, positions=[1, 6, 10])
    )
])
def test_empty(reference, query):
    result = getSut().align(reference, query, __peak(0))

    assert len(result.alignedPairs) == 0


@pytest.mark.parametrize("reference,query,peak", [
    (
            OpticalMap(1, length=9, positions=[0, 4, 8]),
            OpticalMap(1, length=9, positions=[0, 4, 8]),
            __peak(0)
    ),
    (
            OpticalMap(1, length=10, positions=[0, 5, 9]),
            OpticalMap(1, length=10, positions=[0, 5, 9]),
            __peak(0)
    ),
    (
            OpticalMap(1, length=30, positions=[10, 15, 19]),
            OpticalMap(1, length=10, positions=[0, 5, 9]),
            __peak(10)
    )
])
def test_perfectMatch(reference, query, peak):
    result = getSut().align(reference, query, peak)

    assert result.alignedPairs == [(1, 1), (2, 2), (3, 3)]


def test_perfectMatch_reverseStrand():
    reference = OpticalMap(1, length=9, positions=[0, 4, 8])
    query = OpticalMap(1, length=9, positions=[0, 4, 8])

    result = getSut().align(reference, query, __peak(), True)

    assert result.alignedPairs == [(1, 3), (2, 2), (3, 1)]


def test_ignoresExtraPositionsOnReferenceBeforeAndAfterAlignment():
    reference = OpticalMap(1, length=30, positions=[9, 10, 15, 19, 20])
    query = OpticalMap(1, length=10, positions=[0, 5, 9])

    result = getSut().align(reference, query, __peak(10))

    assert result.referenceStartPosition == 10
    assert result.alignedPairs == [(2, 1), (3, 2), (4, 3)]


def test_alignsPositionsWithinMaxDistanceOnly():
    maxDistance = 10
    reference = OpticalMap(1, length=300, positions=[89, 115, 150, 185, 210])
    query = OpticalMap(1, length=100, positions=[0, 25, 50, 75, 99])

    result = getSut(maxDistance).align(reference, query, __peak(100))

    assert result.alignedPairs == [(2, 2), (3, 3), (4, 4)]


def test_returnsCorrectQueryShifts():
    maxDistance = 10
    reference = OpticalMap(1, length=300, positions=[105, 115, 150, 185, 194])
    query = OpticalMap(1, length=100, positions=[0, 25, 50, 75, 99])

    result = getSut(maxDistance).align(reference, query, __peak(100))

    assert list(map(AlignedPair.queryShiftSelector, result.alignedPairs)) == [-5, 10, 0, -10, 5]


def test_ignoresPositionBeyondMaxDistance():
    reference = OpticalMap(1, length=300, positions=[99, 149, 199])
    query = OpticalMap(1, length=100, positions=[0, 49, 89])
    maxDistance = 10

    result = getSut(maxDistance).align(reference, query, __peak(99))

    assert result.referenceEndPosition == 149
    assert result.alignedPairs == [(1, 1), (2, 2)]


@pytest.mark.parametrize("reference,query", [
    (
            OpticalMap(1, length=300, positions=[100, 149, 174, 189]),
            OpticalMap(1, length=100, positions=[0, 49, 89])
    )
])
def test_handlesDeletions_referencePosition3OutOfRange_NotAligned(reference, query):
    result = getSut(10).align(reference, query, __peak(100))

    assert result.referenceEndPosition == 189
    assert result.alignedPairs == [(1, 1), (2, 2), (4, 3)]


@pytest.mark.parametrize("reference,query,expected", [
    (
            OpticalMap(1, length=300, positions=[100, 149, 150, 189]),
            OpticalMap(1, length=100, positions=[0, 49, 89]),
            [(1, 1), (2, 2), (3, None), (4, 3)]
    ), (
            OpticalMap(1, length=300, positions=[100, 149, 188, 189]),
            OpticalMap(1, length=100, positions=[0, 49, 89]),
            [(1, 1), (2, 2), (3, None), (4, 3)]
    )
])
def test_handlesDeletions_atReferencePosition3(reference, query, expected):
    result = getSut(10).align(reference, query, __peak(100))

    assert result.referenceEndPosition == 189
    assert result.positions == expected


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
    result = getSut(10).align(reference, query, __peak(100))

    assert result.referenceEndPosition == 189
    assert result.alignedPairs == [(1, 1), (2, 2), (3, 4)]


@pytest.mark.parametrize("reference,query,refStart,refEnd,queryStart,queryEnd", [
    (
            OpticalMap(1, length=101, positions=[0, 100]),
            OpticalMap(1, length=101, positions=[0, 100]),
            0, 100, 0, 100
    ), (
            OpticalMap(1, length=100, positions=[1, 100]),
            OpticalMap(1, length=100, positions=[1, 100]),
            1, 100, 1, 100
    ), (
            OpticalMap(1, length=200, positions=[1, 105, 200]),
            OpticalMap(1, length=100, positions=[5, 100]),
            1, 105, 5, 100
    )
])
def test_correctAlignmentStartAndEndPositions(reference, query, refStart, refEnd, queryStart, queryEnd):
    result = getSut(10).align(reference, query, __peak(0))

    assert result.referenceStartPosition == refStart
    assert result.referenceEndPosition == refEnd
    assert result.queryStartPosition == queryStart
    assert result.queryEndPosition == queryEnd


@pytest.mark.parametrize("query,positions", [
    (
            OpticalMap(1, length=12, positions=[3, 5, 11]),
            [(1, None), (None, 1), (2, 2), (3, None), (None, 3)]
    ), (
            OpticalMap(1, length=12, positions=[1, 3, 5]),
            [(1, 1), (None, 2), (2, 3), (3, None)]
    ), (
            OpticalMap(1, length=12, positions=[1, 9, 10]),
            [(1, 1), (2, None), (3, 2), (None, 3)]
    ), (
            OpticalMap(1, length=12, positions=[1, 4, 6, 10]),
            [(1, 1), (None, 2), (2, None), (None, 3), (3, None), (None, 4)]
    )
])
def test_returnsUnmatchedPositions(query, positions: List[AlignmentPosition]):
    reference = OpticalMap(1, length=10, positions=[1, 5, 9])
    result = getSut().align(reference, query, __peak(0))

    assert result.positions == positions


def test_multiplePeaks():
    reference = OpticalMap(1, length=111, positions=[0, 4, 100, 110])
    query = OpticalMap(1, length=61, positions=[0, 4, 50, 60])
    peaks = [__peak(0), __peak(50)]
    result = getSut().align(reference, query, peaks)

    assert len(result.segments) == 2
    assert result.segments[0].positions == [(1, 1), (2, 2), (None, 3), (None, 4)]
    assert result.segments[1].positions == [(None, 1), (None, 2), (3, 3), (4, 4)]


if __name__ == '__main__':
    pytest.main(args=[__file__])
