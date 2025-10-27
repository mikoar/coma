from typing import List

import pytest

from src.alignment.segment_chainer import SegmentChainer, SequentialityScorer
from src.alignment.segments import AlignmentSegment
from tests.test_doubles.alignment_segment_builder import AlignmentSegmentBuilder
from tests.test_doubles.scored_aligned_pair_builder import ScoredAlignedPairBuilder
import itertools

chainer = SegmentChainer(SequentialityScorer(1., 0))


@pytest.mark.parametrize("segments, expectedOrder", [
    pytest.param([
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(400).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(500).withQueryPosition(600).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(700).withQueryPosition(700).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(800).withQueryPosition(800).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(100).withQueryPosition(100).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(300).build())
            .withScore(1200.)
            .build()
    ], [2, 3, 1]),
    pytest.param([
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(400).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(500).withQueryPosition(600).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(200).withQueryPosition(300).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(400).withQueryPosition(300).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(100).withQueryPosition(100).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(300).build())
            .withScore(1200.)
            .build()
    ], [3, 2, 1]),
])
def test_chain(segments: List[AlignmentSegment], expectedOrder: List[int]):
    qPositions = list(itertools.chain(*[s.getQueryLabels().positions for s in segments]))
    chainedSegments = chainer.chain(segments, sorted(qPositions))
    assert chainedSegments == [s for _, s in sorted(zip(expectedOrder, segments), key=lambda x: x[0])]


def test_chain_dropsOutlyingSegment():
    segments = [
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(100).withQueryPosition(100).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(300).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(400).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(500).withQueryPosition(600).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(100).withQueryPosition(300).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(500).withQueryPosition(600).build())
            .withScore(1200.)
            .build()
    ]
    qPositions = list(itertools.chain(*[s.getQueryLabels().positions for s in segments]))
    chainedSegments = chainer.chain(segments, sorted(qPositions))
    assert chainedSegments == [segments[0], segments[1]]


def test_chain_preservesEmptySegments():
    segments = [
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(100).withQueryPosition(100).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(300).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(300).withQueryPosition(400).build())
            .withPosition(ScoredAlignedPairBuilder().withReferencePosition(500).withQueryPosition(600).build())
            .withScore(1200.)
            .build(),
        AlignmentSegmentBuilder()
            .withScore(0.)
            .build()
    ]
    qPositions = list(itertools.chain(*[s.getQueryLabels().positions for s in segments]))
    chainedSegments = chainer.chain(segments, sorted(qPositions))
    assert len(chainedSegments) == 3


def test_chain_emptySegmentsOnly():
    segments = [
        AlignmentSegmentBuilder()
            .withScore(0.)
            .build()
    ]
    qPositions = list(itertools.chain(*[s.getQueryLabels().positions for s in segments]))
    chainedSegments = chainer.chain(segments, sorted(qPositions))
    assert chainedSegments == segments


if __name__ == '__main__':
    pytest.main(args=[__file__])
