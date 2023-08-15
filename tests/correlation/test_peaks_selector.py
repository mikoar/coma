import pytest

from src.correlation.optical_map import EmptyInitialAlignment, OpticalMap
from src.correlation.peaks_selector import PeaksSelector


def test_noPeaks_returnsEmpty():
    correlations = [EmptyInitialAlignment(OpticalMap(1, 1, []), OpticalMap(1, 1, []), 0, 0)]
    peaks = PeaksSelector(1).selectPeaks(correlations)
    assert len(peaks) == 0


if __name__ == '__main__':
    pytest.main(args=[__file__])
