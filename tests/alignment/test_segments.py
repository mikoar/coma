from __future__ import annotations

from typing import List, Tuple

import pytest

from src.alignment.segments import AlignmentSegment, AlignmentSegmentsWithResolvedConflicts
from tests.test_doubles.alignment_segment_stub import ScoredAlignedPairStub, \
    ScoredNotAlignedPositionStub


def __segment(scoredPositionTuples: List[Tuple[int | None, int | None, float]]):
    positions = list(map(
        lambda pair: ScoredNotAlignedPositionStub(pair[0], pair[1], pair[2]) if not pair[0] or not pair[1] else
        ScoredAlignedPairStub(pair[0], pair[1], 0, pair[2]), scoredPositionTuples))
    return AlignmentSegment(positions, sum(p.score for p in positions))


@pytest.mark.parametrize("segment0, segment1", [
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
        __segment([(4, 4, 100.), (5, 5, 100.), (6, 6, 100.)]),
        id="segment0 preceding segment 1"
    ),
    pytest.param(
        __segment([(4, 4, 100.), (5, 5, 100.), (6, 6, 100.)]),
        __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
        id="segment0 following segment 1"
    )
])
def test_conflicts_noConflicts_returnsUnchanged(segment0, segment1):
    segmentsAfter = AlignmentSegmentsWithResolvedConflicts.create([segment0, segment1]).segments

    assert segmentsAfter[0] == segment0
    assert segmentsAfter[1] == segment1


@pytest.mark.parametrize("segment0, segment1, expectedSegments", [
    (
            __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
            AlignmentSegment.empty,
            [__segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])]
    ), (
            AlignmentSegment.empty,
            __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
            [__segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])]
    ),
])
def test_conflicts_withEmptySegment_returnsUnchanged(
        segment0, segment1, expectedSegments: List[AlignmentSegment]):
    segmentsAfter = AlignmentSegmentsWithResolvedConflicts.create([segment0, segment1]).segments

    assert segmentsAfter == expectedSegments


checkForConflicts_overlappingSegments_trimsWhileMaxingScore_params = [
    pytest.param(
        __segment([(1, 1, 101.)]),
        __segment([(1, 2, 100.)]),
        [__segment([(1, 1, 101.)])],
        id="single position segments ref conflict"
    ),
    pytest.param(
        __segment([(1, 2, 100.), (2, 3, 100.), (3, 4, 100.)]),
        __segment([(1, 1, 101.)]),
        [__segment([(2, 3, 100.), (3, 4, 100.)]),
         __segment([(1, 1, 101.)])],
        id="single position segment and a larger segment ref conflict at the start, larger segment gets trimmed"
    ),
    pytest.param(
        __segment([(1, 2, 101.), (2, 3, 100.), (3, 4, 100.)]),
        __segment([(1, 1, 100.)]),
        [__segment([(1, 2, 101.), (2, 3, 100.), (3, 4, 100.)])],
        id="single position segment and a larger segment ref conflict at the start, "
           "single position segment gets dropped"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
        __segment([(3, 4, 101.)]),
        [__segment([(1, 1, 100.), (2, 2, 100.)]),
         __segment([(3, 4, 101.)])],
        id="single position segment and a larger segment ref conflict at the end, larger segment gets trimmed"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 101.)]),
        __segment([(3, 4, 100.)]),
        [__segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 101.)])],
        id="single position segment and a larger segment ref conflict at the end, single position segment gets dropped"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.)]),
        __segment([(2, 2, 101.), (3, 3, 100.)]),
        [__segment([(1, 1, 100.)]),
         __segment([(2, 2, 101.), (3, 3, 100.)])],
        id="ref and query conflict at (2,2)"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.)]),
        __segment([(2, 3, 101.), (3, 4, 100.)]),
        [__segment([(1, 1, 100.)]),
         __segment([(2, 3, 101.), (3, 4, 100.)])],
        id="ref conflict at (2,2) and (2,3)"
    ),
    pytest.param(
        __segment([(1, 1, 100.), (2, 2, 100.)]),
        __segment([(3, 2, 101.), (4, 3, 100.)]),
        [__segment([(1, 1, 100.)]),
         __segment([(3, 2, 101.), (4, 3, 100.)])],
        marks=pytest.mark.xfail(reason="not implemented yet"),
        id="query conflict at (2,2) and (3,2)"
    ),
]


@pytest.mark.parametrize("segment0, segment1, expectedSegments",
                         checkForConflicts_overlappingSegments_trimsWhileMaxingScore_params)
def test_conflicts_overlappingSegments_trimsWhileMaxingScore(
        segment0, segment1, expectedSegments: List[AlignmentSegment]):
    segmentsAfter = AlignmentSegmentsWithResolvedConflicts.create([segment0, segment1]).segments
    assert segmentsAfter == expectedSegments


@pytest.mark.parametrize("segment0, segment1, expectedSegments",
                         checkForConflicts_overlappingSegments_trimsWhileMaxingScore_params)
def test_conflicts_overlappingSegments_trimsWhileMaxingScore_reversed(
        segment0, segment1, expectedSegments: List[AlignmentSegment]):
    segmentsAfter = AlignmentSegmentsWithResolvedConflicts.create([segment1, segment0]).segments
    assert segmentsAfter == expectedSegments[::-1]


def test_sliceByReference():
    segment = __segment(
        [(10, 1, 100.), (11, None, -50.), (None, 2, -50.), (12, 3, 100.), (None, 4, -50.), (13, 5, 100.),
         (None, 6, -50.), (14, 7, 100.), (15, 8, 100.)])
    sliced = segment.sliceByReference(12, 14)
    assert sliced.positions == [(12, 3), (None, 4), (13, 5)]
    assert sliced.segmentScore == 150.


def test_subtract():
    segment0 = __segment([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.), (4, 4, 100.)])
    segment1 = __segment([(2, 2, 100.), (3, 3, 100.)])
    result = segment0 - segment1
    assert result == __segment([(1, 1, 100.), (4, 4, 100.)])


def test_subtract_noPositionsLeft_returnsEmptySegment():
    segment0 = __segment([(1, 1, 100.), (2, 2, 100.)])
    segment1 = __segment([(1, 1, 100.), (2, 2, 100.)])
    result = segment0 - segment1
    assert result == AlignmentSegment.empty


if __name__ == '__main__':
    pytest.main(args=[__file__])
