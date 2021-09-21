
import pytest

from src.alignment_score import SegmentScorer


def test_fullAlignment():
    alignmentScores = [5., 6., 7., 4., 3.]
    segment = SegmentScorer().getSegmentWithMaxScore(alignmentScores)
    assert 25. == segment.segmentScore
    assert 0 == segment.start
    assert 4 == segment.end


def test_partialAlignment():
    alignmentScores = [-1., 5., 6., 7., 4., 3., -1]
    segment = SegmentScorer().getSegmentWithMaxScore(alignmentScores)
    assert 25. == segment.segmentScore
    assert 1 == segment.start
    assert 5 == segment.end


def test_alignmentWithGap():
    alignmentScores = [1., 1., -3., 2., 1., -3., 2.]
    segment = SegmentScorer().getSegmentWithMaxScore(alignmentScores)
    assert 3. == segment.segmentScore
    assert 3 == segment.start
    assert 4 == segment.end


if __name__ == '__main__':
    pytest.main(args=[__file__])
