from __future__ import annotations

from typing import List, Tuple

import pytest

from src.alignment.alignment_position import ScoredAlignedPair
from src.alignment.segments import AlignmentSegmentsFactory, AlignmentSegment, AlignmentSegmentsWithoutConflicts
from tests.test_doubles.alignment_segment_stub import AlignedPairStub, ScoredAlignedPairStub, \
    ScoredNotAlignedPositionStub


def __segment(positionTuples: List[Tuple[int | None, int | None, float]]):
    positions = list(map(
        lambda pair: ScoredNotAlignedPositionStub(pair[0], pair[1], pair[2]) if not pair[0] or not pair[1] else
        ScoredAlignedPairStub(pair[0], pair[1], 0, pair[2]), positionTuples))
    return AlignmentSegment(positions, sum(p.score for p in positions))


def __scoredPositions(scores: List[float]):
    return list(map(lambda score: ScoredAlignedPair(AlignedPairStub(0, 0), score), scores))


def test_getSegments_fullAlignment():
    positions = __scoredPositions([5., 6., 7., 4., 3.])
    segments = AlignmentSegmentsFactory(5, 10).getSegments(positions)
    assert len(segments) == 1
    segment = segments[0]
    assert segment.segmentScore == 25.
    assert segment.positions[0].score == positions[0].score
    assert segment.positions[-1].score == positions[4].score
    assert len(segment.positions) == 5


def test_getSegments_empty():
    segments = AlignmentSegmentsFactory(5, 10).getSegments([])
    assert len(segments) == 0


def test_getSegments_partialAlignment():
    positions = __scoredPositions([-1., 5., 6., 7., 4., 3., -1])
    segment = AlignmentSegmentsFactory(5, 10).getSegments(positions)[0]
    assert segment.segmentScore == 25.
    assert segment.positions[0].score == positions[1].score
    assert segment.positions[-1].score == positions[5].score


def test_getSegments_alignmentWithGap():
    positions = __scoredPositions([1., 1., -3., 2., 1., -3., 2.])
    segment = AlignmentSegmentsFactory(3, 1).getSegments(positions)[0]
    assert segment.segmentScore == 3.
    assert segment.positions[0].score == positions[3].score
    assert segment.positions[-1].score == positions[4].score


def test_getSegments_multipleSegments():
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


def test_getSegments_filterSegment():
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


def test_AlignmentSegments_noConflicts_returnsUnchanged():
    segmentBefore0 = __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])
    segmentBefore1 = __segment([(4, 4, 100.), (5, 5, 100.), (6, 6, 100.)])

    segmentsAfter = AlignmentSegmentsWithoutConflicts.create([segmentBefore0, segmentBefore1]).segments

    assert segmentsAfter[0] == segmentBefore0
    assert segmentsAfter[1] == segmentBefore1


@pytest.mark.parametrize("segment0, segment1, expectedSegment0, expectedSegment1", [
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.)]),
        __segment([(2, 2, 101.), (3, 3, 100.)]),
        __segment([(1, 1, 100.)]),
        __segment([(2, 2, 101.), (3, 3, 100.)]),
        id="conflict at (2,2)"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.)]),
        __segment([(2, 3, 101.), (3, 4, 100.)]),
        __segment([(1, 1, 100.)]),
        __segment([(2, 3, 101.), (3, 4, 100.)]),
        id="ref conflict at (2,2) and (2,3)"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.)]),
        __segment([(3, 2, 101.), (4, 3, 100.)]),
        __segment([(1, 1, 100.)]),
        __segment([(3, 2, 101.), (4, 3, 100.)]),
        id="query conflict at (2,2) and (3,2)"
    ),
])
def test_AlignmentSegments_overlappingSegments_trimsWhileMaxingScore(
        segment0, segment1, expectedSegment0, expectedSegment1):
    segmentsAfter = AlignmentSegmentsWithoutConflicts.create([segment0, segment1]).segments

    assert segmentsAfter[0] == expectedSegment0
    assert segmentsAfter[1] == expectedSegment1


def test_sliceByReference():
    segment = __segment(
        [(10, 1, 100.), (11, None, -50.), (None, 2, -50.), (12, 3, 100.), (None, 4, -50.), (13, 5, 100.),
         (None, 6, -50.), (14, 7, 100.), (15, 8, 100.)])
    sliced = segment.sliceByReference(12, 14)
    assert sliced.positions == [(12, 3), (None, 4), (13, 5)]
    assert sliced.segmentScore == 150.


if __name__ == '__main__':
    pytest.main(args=[__file__])
