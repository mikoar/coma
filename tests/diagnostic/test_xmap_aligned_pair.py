import pytest

from src.diagnostic.xmap_alignment import XmapAlignedPair, XmapAlignmentPosition, XmapAlignedPairWithDistance


def __getPair(referencePosition: int, queryPosition: int):
    return XmapAlignedPair(XmapAlignmentPosition(0, referencePosition), XmapAlignmentPosition(0, queryPosition))


def test_calculateDistance_firstPair_distanceShouldBe0():
    pairWithDistance = XmapAlignedPairWithDistance.calculateDistance(__getPair(123, 456), None, False)
    assert pairWithDistance.distance == 0


@pytest.mark.parametrize("pair, expectedDistance", [
    (__getPair(200, 150), 0),
    (__getPair(210, 150), -10),
    (__getPair(190, 150), 10)
])
def test_calculateDistance_distanceIsRelativeToFirstPair(pair: XmapAlignedPair, expectedDistance: int):
    pairWithDistance = XmapAlignedPairWithDistance.calculateDistance(pair, __getPair(100, 50), False)
    assert pairWithDistance.distance == expectedDistance


@pytest.mark.parametrize("pair, expectedDistance", [
    (__getPair(110, 40), 0),
    (__getPair(120, 40), -10),
    (__getPair(100, 40), 10),
])
def test_calculateDistance_distanceIsRelativeToFirstPair_reverseStrand(pair: XmapAlignedPair, expectedDistance: int):
    pairWithDistance = XmapAlignedPairWithDistance.calculateDistance(pair, __getPair(100, 50), True)
    assert pairWithDistance.distance == expectedDistance


if __name__ == '__main__':
    pytest.main(args=[__file__])
