import pytest

from src.diagnostic.benchmark_alignment import BenchmarkAlignedPair, BenchmarkAlignmentPosition, \
    BenchmarkAlignedPairWithDistance


def __getPair(referencePosition: int, queryPosition: int):
    return BenchmarkAlignedPair(BenchmarkAlignmentPosition(0, referencePosition), BenchmarkAlignmentPosition(0, queryPosition))


def test_calculateDistance_firstPair_distanceShouldBe0():
    pairWithDistance = BenchmarkAlignedPairWithDistance.calculateDistance(__getPair(123, 456), None, False)
    assert pairWithDistance.distance == 0


@pytest.mark.parametrize("pair, expectedDistance", [
    (__getPair(200, 150), 0),
    (__getPair(210, 150), -10),
    (__getPair(190, 150), 10)
])
def test_calculateDistance_distanceIsRelativeToFirstPair(pair: BenchmarkAlignedPair, expectedDistance: int):
    pairWithDistance = BenchmarkAlignedPairWithDistance.calculateDistance(pair, __getPair(100, 50), False)
    assert pairWithDistance.distance == expectedDistance


@pytest.mark.parametrize("pair, expectedDistance", [
    (__getPair(110, 40), 0),
    (__getPair(120, 40), -10),
    (__getPair(100, 40), 10),
])
def test_calculateDistance_distanceIsRelativeToFirstPair_reverseStrand(pair: BenchmarkAlignedPair, expectedDistance: int):
    pairWithDistance = BenchmarkAlignedPairWithDistance.calculateDistance(pair, __getPair(100, 50), True)
    assert pairWithDistance.distance == expectedDistance


if __name__ == '__main__':
    pytest.main(args=[__file__])
