from __future__ import annotations

from typing import List

import pytest

from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments import AlignmentSegment, EmptyAlignmentSegment
from tests.test_doubles.alignment_segment_stub import AlignedPairStub, AlignmentSegmentStub
from tests.test_doubles.mock_segment_chainer import MockSegmentChainer


def test_conflicts_noConflicts_returnsUnchanged():
    segment0 = AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])
    segment1 = AlignmentSegmentStub.createFromPairs([(4, 4, 100.), (5, 5, 100.), (6, 6, 100.)])
    resolver = AlignmentSegmentConflictResolver(MockSegmentChainer())

    qPositions = segment0.getQueryLabels().positions + segment1.getQueryLabels().positions 
    segmentsAfter = resolver.resolveConflicts([segment0, segment1], sorted(qPositions)).segments

    assert segmentsAfter[0] == segment0
    assert segmentsAfter[1] == segment1


@pytest.mark.parametrize("segment0, segment1, expectedSegments", [
    (
            AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
            EmptyAlignmentSegment(),
            [AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
             EmptyAlignmentSegment()]
    ), (
            EmptyAlignmentSegment(),
            AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
            [EmptyAlignmentSegment(),
             AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])]
    ),
])
def test_conflicts_withEmptySegment_returnsUnchanged(
        segment0, segment1, expectedSegments: List[AlignmentSegment]):
    resolver = AlignmentSegmentConflictResolver(MockSegmentChainer())
    qPositions = segment0.getQueryLabels().positions + segment1.getQueryLabels().positions 
    segmentsAfter = resolver.resolveConflicts([segment0, segment1], sorted(qPositions)).segments
    assert segmentsAfter == expectedSegments


checkForConflicts_overlappingSegments_trimsWhileMaxingScore_params = [
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 101.)]),
        AlignmentSegmentStub.createFromPairs([(1, 2, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 101.)]),
         EmptyAlignmentSegment()],
        id="single position segments ref conflict"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 101.)]),
        AlignmentSegmentStub.createFromPairs([(1, 2, 100.), (2, 3, 100.), (3, 4, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 101.)]),
         AlignmentSegmentStub.createFromPairs([(2, 3, 100.), (3, 4, 100.)])],
        id="single position segment and a larger segment ref conflict at the start, larger segment gets trimmed"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.)]),
        AlignmentSegmentStub.createFromPairs([(1, 2, 101.), (2, 3, 100.), (3, 4, 100.)]),
        [EmptyAlignmentSegment(),
         AlignmentSegmentStub.createFromPairs([(1, 2, 101.), (2, 3, 100.), (3, 4, 100.)])],
        id="single position segment and a larger segment ref conflict at the start, "
           "single position segment gets dropped"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
        AlignmentSegmentStub.createFromPairs([(3, 4, 101.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.)]),
         AlignmentSegmentStub.createFromPairs([(3, 4, 101.)])],
        id="single position segment and a larger segment ref conflict at the end, larger segment gets trimmed"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 101.)]),
        AlignmentSegmentStub.createFromPairs([(3, 4, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 101.)]),
         EmptyAlignmentSegment()],
        id="single position segment and a larger segment ref conflict at the end, single position segment gets dropped"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.)]),
        AlignmentSegmentStub.createFromPairs([(2, 2, 101.), (3, 3, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 100.)]),
         AlignmentSegmentStub.createFromPairs([(2, 2, 101.), (3, 3, 100.)])],
        id="ref and query conflict at (2,2)"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.)]),
        AlignmentSegmentStub.createFromPairs([(2, 3, 101.), (3, 4, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 100.)]),
         AlignmentSegmentStub.createFromPairs([(2, 3, 101.), (3, 4, 100.)])],
        id="ref conflict at (2,2) and (2,3)"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.)]),
        AlignmentSegmentStub.createFromPairs([(3, 2, 101.), (4, 3, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(1, 1, 100.)]),
         AlignmentSegmentStub.createFromPairs([(3, 2, 101.), (4, 3, 100.)])],
        id="query conflict at (2,2) and (3,2)"
    ),
    pytest.param(
        AlignmentSegmentStub.createFromPairs([(4, 4, 100.), (5, 7, 100.)]),
        AlignmentSegmentStub.createFromPairs([(6, 6, 101.), (8, 8, 100.)]),
        [AlignmentSegmentStub.createFromPairs([(4, 4, 100.)]),
         AlignmentSegmentStub.createFromPairs([(6, 6, 101.), (8, 8, 100.)])],
        id="cross conflict"
    ),
]


@pytest.mark.parametrize("segment0, segment1, expectedSegments",
                         checkForConflicts_overlappingSegments_trimsWhileMaxingScore_params)
def test_conflicts_overlappingSegments_trimsWhileMaxingScore(
        segment0: AlignmentSegment, segment1: AlignmentSegment, expectedSegments: List[AlignmentSegment]):
    resolver = AlignmentSegmentConflictResolver(MockSegmentChainer())
    qPositions = segment0.getQueryLabels().positions + segment1.getQueryLabels().positions 
    segmentsAfter = resolver.resolveConflicts([segment0, segment1], sorted(qPositions)).segments
    assert segmentsAfter == expectedSegments


@pytest.mark.parametrize("start, stop", [
    pytest.param(AlignedPairStub(12, 3), AlignedPairStub(13, 5)),
])
def test_slice(start, stop):
    segment = AlignmentSegmentStub.createFromPairs(
        [(10, 1, 100.), (11, None, -50.), (None, 2, -50.), (12, 3, 100.), (None, 4, -50.), (13, 5, 100.),
         (None, 6, -50.), (14, 7, 100.), (15, 8, 100.)])

    sliced = segment.slice(start, stop)

    assert sliced.positions == [(12, 3), (None, 4), (13, 5)]
    assert sliced.segmentScore == 150.
    assert len(sliced.allPeakPositions) == 9


@pytest.mark.parametrize("start, stop", [
    pytest.param(AlignedPairStub(12, 3), AlignedPairStub(13, 5)),
])
def test_slice_parameterEnd_keepsUnpaired(start, stop):
    segment = AlignmentSegmentStub.createFromPairs(
        [(10, 1, 100.), (11, None, -50.), (None, 2, -50.), (12, 3, 100.), (13, 4, -50.), (None, 5, 100.),
         (None, 6, -50.), (None, 7, 100.), (None, 8, 100.)])

    sliced = segment.slice(start, stop)

    assert sliced.positions == [(12, 3), (13, 4), (None, 5)]
    assert sliced.segmentScore == 150.
    assert len(sliced.allPeakPositions) == 9


def test_subtract():
    segment0 = AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.), (4, 4, 100.)])
    segment1 = AlignmentSegmentStub.createFromPairs([(2, 2, 100.), (3, 3, 100.)])

    result = segment0 - segment1

    assert result == AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (4, 4, 100.)])
    assert len(result.allPeakPositions) == 4


def test_subtract_noPositionsLeft_returnsEmptySegment():
    segment0 = AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.)])
    segment1 = AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.)])

    result = segment0 - segment1

    assert result == EmptyAlignmentSegment()


if __name__ == '__main__':
    pytest.main(args=[__file__])
