from dataclasses import dataclass
from typing import List, Iterator

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
        sortedPeaks = sorted((SelectedPeak(c, p) for c in correlations for p in c.peaks), key=lambda sp: sp.peak.score, reverse=True)
        return sortedPeaks[:self.count]
