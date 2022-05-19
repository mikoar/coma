from __future__ import annotations

from typing import List

import pytest

from src.alignment.alignment_position import ScoredAlignedPair
from src.alignment.segments_factory import AlignmentSegmentsFactory
from tests.test_doubles.alignment_segment_stub import AlignedPairStub, ScoredAlignedPairStub


def __scoredPositions(scores: List[float]):
    return list(map(lambda score: ScoredAlignedPair(AlignedPairStub(0, 0), score), scores))


def test_fullAlignment():
    positions = __scoredPositions([5., 6., 7., 4., 3.])
    segments = AlignmentSegmentsFactory(5, 10).getSegments(positions)
    assert len(segments) == 1
    segment = segments[0]
    assert segment.segmentScore == 25.
    assert segment.positions[0].score == positions[0].score
    assert segment.positions[-1].score == positions[4].score
    assert len(segment.positions) == 5


def test_empty():
    segments = AlignmentSegmentsFactory(5, 10).getSegments([])
    assert len(segments) == 0


def test_partialAlignment():
    positions = __scoredPositions([-1., 5., 6., 7., 4., 3., -1])
    segment = AlignmentSegmentsFactory(5, 10).getSegments(positions)[0]
    assert segment.segmentScore == 25.
    assert segment.positions[0].score == positions[1].score
    assert segment.positions[-1].score == positions[5].score


def test_alignmentWithGap():
    positions = __scoredPositions([1., 1., -3., 2., 1., -3., 2.])
    segment = AlignmentSegmentsFactory(3, 1).getSegments(positions)[0]
    assert segment.segmentScore == 3.
    assert segment.positions[0].score == positions[3].score
    assert segment.positions[-1].score == positions[4].score


def test_multipleSegments():
    positions = __scoredPositions([1., 1., 1., 1., 1., -1., -1., -1., 1., 1.])
    minScore = 2
    threshold = 3
    segments = AlignmentSegmentsFactory(minScore, threshold).getSegments(positions)

    assert len(segments) == 2
    assert segments[0].segmentScore == 5.
    assert segments[0].positions[0].score == positions[0].score
    assert segments[0].positions[-1].score == positions[4].score
    assert segments[1].segmentScore == 2.
    assert segments[0].positions[0].score == positions[8].score
    assert segments[0].positions[-1].score == positions[9].score


def test_filterSegment():
    pairs = list(map(lambda pair: ScoredAlignedPairStub(pair[0], pair[1], 0, pair[2]),
                     [(1616, 1, 122.7), (1617, 2, 25.9), (1618, 3, 97.7), (1619, 4, 197.7), (1621, 5, 166.0),
                      (1622, 6, 155.3), (1623, 7, 159.4), (1624, 8, 94.9), (1625, 9, 165.1), (1626, 10, 22.4),
                      (1627, 11, 34.0), (1630, 14, -562.9), (1632, 16, -1533.2), (1633, 17, -1399.0),
                      (1636, 20, -1400.0), (1637, 21, -872.8), (1639, 22, -1820.0)]))

    segments = AlignmentSegmentsFactory(200, 0).getSegments(pairs)

    assert len(segments) == 1
    assert segments[0].positions == [(1616, 1), (1617, 2), (1618, 3), (1619, 4),
                                     (1621, 5), (1622, 6), (1623, 7), (1624, 8),
                                     (1625, 9), (1626, 10), (1627, 11)]


if __name__ == '__main__':
    pytest.main(args=[__file__])