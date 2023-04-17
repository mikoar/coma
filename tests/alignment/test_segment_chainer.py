from typing import List

import pytest

from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segments import AlignmentSegment
from tests.test_doubles.alignment_segment_stub import AlignmentSegmentStub


@pytest.mark.parametrize("segments, expectedOrder", [
    pytest.param([
        AlignmentSegmentStub.createFromPairs([(3, 4, 80.), (4, 5, 100.), (5, 6, 100.)]),
        AlignmentSegmentStub.createFromPairs([(7, 7, 100.), (8, 8, 100.)]),
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])
    ], [2, 3, 1]),
    pytest.param([
        AlignmentSegmentStub.createFromPairs([(3, (5, 300), 80.), (4, (4, 400), 100.), (5, (3, 500), 100.)]),
        AlignmentSegmentStub.createFromPairs([(7, (2, 600), 100.), (8, (1, 700), 100.)]),
        AlignmentSegmentStub.createFromPairs([(1, (8, 0), 100.), (2, (7, 100), 100.), (3, (6, 200), 100.)])
    ], [2, 3, 1], id="reverse strand"),
    pytest.param([
        AlignmentSegmentStub.createFromPairs([(3, 4, 80.), (4, 5, 100.), (5, 6, 100.)]),
        AlignmentSegmentStub.createFromPairs([(2, 3, 100.), (4, 3, 100.)]),
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)])
    ], [3, 2, 1]),
    pytest.param([
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (None, 3, 0.), (None, 4, 0.)]),
        AlignmentSegmentStub.createFromPairs([(None, 1, 0.), (None, 2, 0.), (3, 3, 100.), (4, 4, 100.)]),
    ], [1, 2]),
])
def test_chain(segments: List[AlignmentSegment], expectedOrder: List[int]):
    chainedSegments = SegmentChainer().chain(segments)
    assert chainedSegments == [s for _, s in sorted(zip(expectedOrder, segments), key=lambda x: x[0])]


def test_chain_dropsOutlyingSegment():
    segments = [
        AlignmentSegmentStub.createFromPairs([(1, 1, 100.), (2, 2, 100.), (3, 3, 100.)]),
        AlignmentSegmentStub.createFromPairs([(3, 4, 80.), (4, 5, 100.), (5, 6, 100.)]),
        AlignmentSegmentStub.createFromPairs([(1, 3, 10.), (2, 5, 120.), (5, 6, 5.)])
    ]
    chainedSegments = SegmentChainer().chain(segments)
    assert chainedSegments == [segments[0], segments[1]]


if __name__ == '__main__':
    pytest.main(args=[__file__])
