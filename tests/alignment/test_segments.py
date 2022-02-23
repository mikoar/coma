from __future__ import annotations

from typing import List

import pytest

from src.alignment.alignment_position import nullAlignedPair, ScoredAlignedPair
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segments import AlignmentSegments
from tests.alignment.test_alignment_position import TestAlignedPair


def test_fullAlignment():
    positions = __scoredPositions([5., 6., 7., 4., 3.])
    segments = AlignmentSegments.getSegments(positions, 5, 10)
    assert len(segments) == 1
    segment = segments[0]
    assert segment.segmentScore == 25.
    assert segment.start == 0
    assert segment.end == 5
    assert len(segment.positions) == 5


def test_partialAlignment():
    positions = __scoredPositions([-1., 5., 6., 7., 4., 3., -1])
    segment = AlignmentSegments.getSegments(positions, 5, 10)[0]
    assert segment.segmentScore == 25.
    assert segment.start == 1
    assert segment.end == 6


def test_alignmentWithGap():
    positions = __scoredPositions([1., 1., -3., 2., 1., -3., 2.])
    segment = AlignmentSegments.getSegments(positions, 3, 1)[0]
    assert segment.segmentScore == 3.
    assert segment.start == 3
    assert segment.end == 5


def test_multipleSegments():
    positions = __scoredPositions([1., 1., 1., 1., 1., -1., -1., -1., 1., 1.])
    minScore = 2
    threshold = 3
    segments = AlignmentSegments.getSegments(positions, minScore, threshold)

    assert len(segments) == 2
    assert segments[0].segmentScore == 5.
    assert segments[0].start == 0
    assert segments[0].end == 5
    assert segments[1].segmentScore == 2.
    assert segments[1].start == 8
    assert segments[1].end == 10


def test_filterSegments():
    row = AlignmentResultRow(TestAlignedPair.createAlignedPairs(
        [(1616, 1, -77.30), (1617, 2, -174.10), (1618, 3, -102.30), (1619, 4, -2.30), (1621, 5, 34.00),
         (1622, 6, -44.70), (1623, 7, -40.60), (1624, 8, 105.10), (1625, 9, 34.90), (1626, 10, 177.60),
         (1627, 11, 166.00), (1630, 14, 762.90), (1632, 16, -1733.20), (1633, 17, 1599.00), (1636, 20, 1600.00),
         (1637, 21, -1072.80), (1639, 22, 2020.00)]))

    filteredRow = AlignmentSegments.filterSegments(row, 200, 2, -200, 1, 0)

    assert filteredRow.alignedPairs == [(1616, 1, -77.30), (1617, 2, -174.10), (1618, 3, -102.30), (1619, 4, -2.30),
                                     (1621, 5, 34.00), (1622, 6, -44.70), (1623, 7, -40.60), (1624, 8, 105.10),
                                     (1625, 9, 34.90), (1626, 10, 177.60), (1627, 11, 166.00)]


def __scoredPositions(scores: List[float]):
    return list(map(lambda score: ScoredAlignedPair(nullAlignedPair, score), scores))


if __name__ == '__main__':
    pytest.main(args=[__file__])
