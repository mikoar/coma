from typing import List

import pytest

from src.correlation.peak import Peak
from src.correlation.peaks_selector import PeaksSelector
from tests.test_doubles.initial_alignment_builder import InitialAlignmentBuilder


def test_noPeaks_returnsEmpty():
    correlations = [InitialAlignmentBuilder().build(), InitialAlignmentBuilder().build()]
    peaks = PeaksSelector(1).selectPeaks(correlations)
    assert len(peaks) == 0


@pytest.mark.parametrize("desiredNumberOfPeaks", [1, 2])
def test_selectsOnePeakFromSingleCorrelationWithOnePeak(desiredNumberOfPeaks: int):
    peak = Peak(100, 50, score=40)
    correlation = InitialAlignmentBuilder().withPeak(peak).build()

    peaks = PeaksSelector(desiredNumberOfPeaks).selectPeaks([correlation])

    assert len(peaks) == 1
    assert peaks[0].peak == peak
    assert peaks[0].primaryCorrelation == correlation


@pytest.mark.parametrize("numberOfPeaks, expectedScores", [
    (0, []),
    (1, [40]),
    (2, [40, 30]),
    (3, [40, 30, 20]),
    (4, [40, 30, 20, 10])
])
def test_selectsMultiplePeaksFromOneCorrelationWithMultiplePeaks_sortedByScore(
        numberOfPeaks: int, expectedScores: List[float]):
    peak1 = Peak(1, 30, score=20)
    peak2 = Peak(2, 20, score=10)
    peak3 = Peak(3, 40, score=30)
    peak4 = Peak(4, 50, score=40)

    correlation = InitialAlignmentBuilder() \
        .withPeak(peak1) \
        .withPeak(peak2) \
        .withPeak(peak3) \
        .withPeak(peak4) \
        .build()

    peaks = PeaksSelector(numberOfPeaks).selectPeaks([correlation])

    assert len(peaks) == numberOfPeaks
    for peak, expectedScore in zip(peaks, expectedScores):
        assert peak.peak.score == expectedScore


def test_selectsMultiplePeaksFromMultipleCorrelations():
    peak11 = Peak(1, 30, score=20)
    peak12 = Peak(2, 20, score=40)
    peak13 = Peak(3, 40, score=30)
    peak21 = Peak(1, 150, score=10)
    peak22 = Peak(2, 140, score=50)

    correlation1 = InitialAlignmentBuilder() \
        .withPeak(peak11) \
        .withPeak(peak12) \
        .withPeak(peak13) \
        .build()

    correlation2 = InitialAlignmentBuilder() \
        .withPeak(peak21) \
        .withPeak(peak22) \
        .build()

    peaks = PeaksSelector(2).selectPeaks([correlation1, correlation2])

    assert len(peaks) == 2
    assert peaks[0].peak == peak22
    assert peaks[1].peak == peak12


if __name__ == '__main__':
    pytest.main(args=[__file__])
