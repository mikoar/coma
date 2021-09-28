from __future__ import annotations

import pytest

from src.alignment.region_scores import RegionScores


def test_fullAlignment():
    alignmentScores = [5., 6., 7., 4., 3.]
    segment = RegionScores(alignmentScores).getSegmentWithMaxScore()
    assert 25. == segment.score
    assert 0 == segment.start
    assert 4 == segment.end


def test_partialAlignment():
    alignmentScores = [-1., 5., 6., 7., 4., 3., -1]
    segment = RegionScores(alignmentScores).getSegmentWithMaxScore()
    assert 25. == segment.score
    assert 1 == segment.start
    assert 5 == segment.end


def test_alignmentWithGap():
    alignmentScores = [1., 1., -3., 2., 1., -3., 2.]
    segment = RegionScores(alignmentScores).getSegmentWithMaxScore()
    assert 3. == segment.score
    assert 3 == segment.start
    assert 4 == segment.end


if __name__ == '__main__':
    pytest.main(args=[__file__])
