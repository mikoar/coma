
import pytest

from vectorise import blur, vectorisePositions, vectoriseSegments


def test_vectorise_simple():
    seq = [2, 5, 8]
    result = list(vectorisePositions(seq, 1))

    assert result == [0, 0, 1, 0, 0, 1, 0, 0, 1]


def test_vectorise_in_the_middle_of_window():
    seq = [5, 15, 35]
    result = list(vectorisePositions(seq, 10))

    assert result == [1, 1, 0, 1]


def test_vectorise_clustered():
    seq = [10, 11, 12, 13, 14, 21]
    result = list(vectorisePositions(seq, 5))

    assert result == [0, 0, 1, 0, 1]


def test_vectoriseSegments():
    segments = [0, 3, 3, 2, 1]
    result = list(vectoriseSegments(segments, 1))

    assert result == [1, 0, 0, 1, 0, 0, 1, 0, 1, 1]


def test_blur_radius_0():
    vector = [0, 0, 1, 0, 0]
    expect = [0, 0, 1, 0, 0]

    result = blur(vector, 0).tolist()

    assert result == expect


def test_blur_radius_1():
    vector = [1, 0, 0, 0, 1, 0, 1, 0, 0]
    expect = [1, 1, 0, 1, 1, 1, 1, 1, 0]

    result = blur(vector, 1).tolist()

    assert result == expect


def test_blur_radius_2():
    vector = [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1]
    expect = [0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1]

    result = blur(vector, 2).tolist()

    assert result == expect


if __name__ == '__main__':
    pytest.main(args=[__file__])
