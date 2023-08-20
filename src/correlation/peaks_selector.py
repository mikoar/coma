from collections import Iterator
from dataclasses import dataclass
from typing import List

from src.correlation.optical_map import InitialAlignment
from src.correlation.peak import Peak


@dataclass
class SelectedPeak:
    primaryCorrelation: InitialAlignment
    peak: Peak


class PeaksSelector:
    def __init__(self, count: int):
        self.count = count

    def selectPeaks(self, correlations: Iterator[InitialAlignment]) -> List[SelectedPeak]:
        peaks = (SelectedPeak(c, p) for c in correlations for p in c.peaks)
        topPeaks = sorted(peaks, key=lambda sp: sp.peak.score, reverse=True)[0:self.count]
        return topPeaks
