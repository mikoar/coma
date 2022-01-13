from __future__ import annotations

import pytest

from src.alignment.region_scores import RegionScores


def test_fullAlignment():
    alignmentScores = [5., 6., 7., 4., 3.]
    segments = RegionScores(alignmentScores).getSegments(5, 10)
    assert len(segments) == 1
    segment = segments[0]
    assert segment.score == 25.
    assert segment.start == 0
    assert segment.end == 4


def test_partialAlignment():
    alignmentScores = [-1., 5., 6., 7., 4., 3., -1]
    segment = RegionScores(alignmentScores).getSegments(5, 10)[0]
    assert segment.score == 25.
    assert segment.start == 1
    assert segment.end == 5


def test_alignmentWithGap():
    alignmentScores = [1., 1., -3., 2., 1., -3., 2.]
    segment = RegionScores(alignmentScores).getSegments(3, 1)[0]
    assert segment.score == 3.
    assert segment.start == 3
    assert segment.end == 4


def test_multipleSegments():
    alignmentScores = [1., 1., 1., 1., 1., -1., -1., -1., 1., 1.]
    minScore = 2
    threshold = 3
    segments = RegionScores(alignmentScores).getSegments(minScore, threshold)

    assert len(segments) == 2
    assert segments[0].score == 5.
    assert segments[0].start == 0
    assert segments[0].end == 4
    assert segments[1].score == 2.
    assert segments[1].start == 8
    assert segments[1].end == 9


if __name__ == '__main__':
    pytest.main(args=[__file__])
