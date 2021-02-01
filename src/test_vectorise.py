
import pytest
from vectorise import vectorise


def test_vectorise_simple():
    seq = [2, 5, 8]
    result = list(vectorise(seq, 1))

    assert result == [0, 0, 1, 0, 0, 1, 0, 0, 1]


def test_vectorise_in_the_middle_of_window():
    seq = [5, 15, 35]
    result = list(vectorise(seq, 10))

    assert result == [1, 1, 0, 1]


def test_vectorise_clustered():
    seq = [10, 11, 12, 13, 14, 21]
    result = list(vectorise(seq, 5))

    assert result == [0, 0, 1, 0, 1]


if __name__ == '__main__':
    pytest.main(args=[__file__])

# TODO: testy na genomie, wyciąć fragment i dopasować, symulowane dane
