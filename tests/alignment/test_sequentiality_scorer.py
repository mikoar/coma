import math

import pytest

from src.alignment.segment_chainer import SequentialityScorer
from tests.test_doubles.alignment_segment_builder import AlignmentSegmentBuilder
from tests.test_doubles.scored_aligned_pair_builder import ScoredAlignedPairBuilder


@pytest.mark.parametrize("pos1, pos2, pos3, pos4, expectedScore", [
    pytest.param((100, 100), (500, 500), (500, 500), (1000, 1000), 0, id="no distance"),
    pytest.param((100, 100), (500, 500), (501, 501), (1000, 1000), -2, id="1 distance on both"),
    pytest.param((100, 100), (500, 500), (600, 600), (1000, 1000), -200, id="100 distance on both"),
    pytest.param((100, 100), (500, 500), (400, 400), (1000, 1000), -200, id="100 overlap on both"),
    pytest.param((100, 100), (500, 500), (500, 600), (1000, 1000), -200, id="100 distance on query"),
    pytest.param((100, 100), (500, 500), (600, 500), (1000, 1000), -200, id="100 distance on reference"),
    pytest.param((100, 100), (500, 500), (500, 400), (1000, 1000), -200, id="100 overlap on query"),
    pytest.param((100, 100), (500, 500), (400, 500), (1000, 1000), -200, id="100 overlap on reference"),
    pytest.param((100, 100), (500, 500), (400, 600), (1000, 1000), -200, id="100 reference overlap and query distance"),
    pytest.param((100, 100), (500, 500), (600, 400), (1000, 1000), -200, id="100 query overlap and reference distance"),
    pytest.param((100, 100), (500, 500), (700, 700), (1000, 1000), -400, id="200 distance on both"),
    pytest.param((100, 100), (500, 500), (300, 300), (1000, 1000), -400, id="200 overlap on both"),
    pytest.param((100, 100), (500, 500), (500, 700), (1000, 1000), -400, id="200 distance on query"),
    pytest.param((100, 100), (500, 500), (700, 500), (1000, 1000), -400, id="200 distance on reference"),
    pytest.param((100, 100), (500, 500), (500, 300), (1000, 1000), -400, id="200 overlap on query"),
    pytest.param((100, 100), (500, 500), (300, 500), (1000, 1000), -400, id="200 overlap on reference"),
    pytest.param((100, 100), (500, 500), (300, 700), (1000, 1000), -400, id="200 reference overlap and query distance"),
    pytest.param((100, 100), (500, 500), (700, 300), (1000, 1000), -400, id="200 query overlap and reference distance"),
    pytest.param((100, 100), (500, 500), (500, 200), (1000, 1000), -math.inf, id="overlap over threshold on query"),
    pytest.param((100, 100), (500, 500), (200, 500), (1000, 1000), -math.inf, id="overlap over threshold on reference"),
])
def test_getScore(pos1, pos2, pos3, pos4, expectedScore):
    previousSegment = AlignmentSegmentBuilder() \
        .withPosition(ScoredAlignedPairBuilder().withReferencePosition(pos1[0]).withQueryPosition(pos1[1]).build()) \
        .withPosition(ScoredAlignedPairBuilder().withReferencePosition(pos2[0]).withQueryPosition(pos2[1]).build()) \
        .build()
    currentSegment = AlignmentSegmentBuilder() \
        .withPosition(ScoredAlignedPairBuilder().withReferencePosition(pos3[0]).withQueryPosition(pos3[1]).build()) \
        .withPosition(ScoredAlignedPairBuilder().withReferencePosition(pos4[0]).withQueryPosition(pos4[1]).build()) \
        .build()

    score = SequentialityScorer(1., 0).getScore(previousSegment, currentSegment)

    assert score == expectedScore
