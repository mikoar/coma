from __future__ import annotations

from typing import List

import pytest

from src.alignment.alignment_position import ScoredAlignedPair
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.correlation.peak import Peak
from tests.test_doubles.alignment_segment_stub import AlignedPairStub, ScoredAlignedPairStub


def __scoredPositions(scores: List[float]):
    return list(map(lambda x: ScoredAlignedPair(AlignedPairStub(x[0], x[0]), x[1]), enumerate(scores)))


def __peak(position: int = 0):
    return Peak(position, 100)


def test_fullAlignment():
    positions = __scoredPositions([5., 6., 7., 4., 3.])
    segments = AlignmentSegmentsFactory(5, 10).getSegments(positions, __peak())
    assert len(segments) == 1
    segment = segments[0]
    assert segment.segmentScore == 25.
    assert segment.positions[0] == positions[0]
    assert segment.positions[-1] == positions[4]
    assert len(segment.positions) == 5
    assert segment.allPeakPositions == positions


def test_noPositions_returnsEmptySegment():
    segments = AlignmentSegmentsFactory(5, 10).getSegments([], __peak())
    assert len(segments) == 1
    emptySegment = segments[0]
    assert len(emptySegment.positions) == 0
    assert len(emptySegment.allPeakPositions) == 0


def test_noQualifyingSegment_returnsEmptySegment():
    positions = __scoredPositions([5., 6., 7., 4., 3.])
    segments = AlignmentSegmentsFactory(99, 10).getSegments(positions, __peak())
    assert len(segments) == 1
    segment = segments[0]
    assert segment.segmentScore == 0
    assert len(segment.positions) == 0
    assert len(segment.allPeakPositions) == 5


def test_partialAlignment():
    positions = __scoredPositions([-1., 5., 6., 7., 4., 3., -1])
    segment = AlignmentSegmentsFactory(5, 10).getSegments(positions, __peak())[0]
    assert segment.segmentScore == 25.
    assert segment.positions[0] == positions[1]
    assert segment.positions[-1] == positions[5]
    assert segment.allPeakPositions == positions


def test_alignmentWithGap():
    positions = __scoredPositions([1., 1., -3., 2., 1., -3., 2.])
    segment = AlignmentSegmentsFactory(3, 1).getSegments(positions, __peak())[0]
    assert segment.segmentScore == 3.
    assert segment.positions[0] == positions[3]
    assert segment.positions[-1] == positions[4]
    assert segment.allPeakPositions == positions


def test_setsPeakPosition():
    positions = __scoredPositions([1.])
    segment = AlignmentSegmentsFactory(1, 1).getSegments(positions, __peak(123))[0]
    assert segment.peak.position == 123


def test_multipleSegments():
    positions = __scoredPositions([1., 1., 1., 1., 1., -1., -1., -1., -1., -1., -1., -1., 1., 1.])
    minScore = 2
    threshold = 3
    segments = AlignmentSegmentsFactory(minScore, threshold).getSegments(positions, __peak())

    assert len(segments) == 2
    assert segments[0].segmentScore == 5.
    assert segments[0].positions[0] == positions[0]
    assert segments[0].positions[-1] == positions[4]
    assert segments[0].allPeakPositions == positions
    assert segments[1].segmentScore == 2.
    assert segments[1].positions[0] == positions[12]
    assert segments[1].positions[-1] == positions[13]
    assert segments[1].allPeakPositions == positions


def test_filterSegment():
    pairs = list(map(lambda pair: ScoredAlignedPairStub(pair[0], pair[1], 0, pair[2]),
                     [(1616, 1, 122.7), (1617, 2, 25.9), (1618, 3, 97.7), (1619, 4, 197.7), (1621, 5, 166.0),
                      (1622, 6, 155.3), (1623, 7, 159.4), (1624, 8, 94.9), (1625, 9, 165.1), (1626, 10, 22.4),
                      (1627, 11, 34.0), (1630, 14, -562.9), (1632, 16, -1533.2), (1633, 17, -1399.0),
                      (1636, 20, -1400.0), (1637, 21, -872.8), (1639, 22, -1820.0)]))

    segments = AlignmentSegmentsFactory(200, 0).getSegments(pairs, __peak())

    assert len(segments) == 1
    assert segments[0].positions == [(1616, 1), (1617, 2), (1618, 3), (1619, 4),
                                     (1621, 5), (1622, 6), (1623, 7), (1624, 8),
                                     (1625, 9), (1626, 10), (1627, 11)]


if __name__ == '__main__':
    pytest.main(args=[__file__])
