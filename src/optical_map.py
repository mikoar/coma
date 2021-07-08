from __future__ import annotations
from functools import cached_property
from scipy.signal import find_peaks
from typing import List
import numpy as np

from typing import Iterable
from scipy.signal import correlate


class OpticalMap:
    def __init__(self, moleculeId: int, sequence: List[int], positions:  List[int], resolution: int) -> None:
        self.sequence = sequence
        self.positions = positions
        self.moleculeId = moleculeId
        self.resolution = resolution

    def correlate(self, reference: OpticalMap, reverseStrand=False, flatten=True):
        if reverseStrand:
            self.sequence.reverse()

        correlation = self.__getCorrelation(reference.sequence, self.sequence)

        if flatten:
            normalizingFactor = self.__getCorrelation(reference.sequence, [1] * len(self.sequence)) + sum(self.sequence)
            correlation = correlation / normalizingFactor

        correlation /= np.max(correlation)

        return CorrelationResult(correlation, self, reference)

    def __getCorrelation(self, reference: Iterable[int], query: Iterable[int]):
        return correlate(reference, query, mode='same')


class CorrelationResult:
    def __init__(self, correlation: List[float], query: OpticalMap, reference: OpticalMap) -> None:
        self.correlation = correlation
        self.query = query
        self.reference = reference

    @cached_property
    def peaks(self):
        return Peaks(self)


class Peaks:
    def __init__(self, correlationResult: CorrelationResult) -> None:
        self.correlationResult = correlationResult
        self.peaks, self.peakProperties = find_peaks(correlationResult.correlation,
                                                     height=0.05,
                                                     prominence=0.2,
                                                     distance=((5 * 10 ** 6) / correlationResult.query.resolution))
        self.score = self.__getScore()
        self.reverseScore = 1 / self.score if self.score else 1
        self.max = self.__getMax()

    def __getScore(self):
        heights = self.__peakHeights

        if len(heights) == 1:
            return 1

        if len(heights) < 2:
            return 0

        twoHighest = sorted(heights, reverse=True)[:2]
        return twoHighest[0] - twoHighest[1]

    def __getMax(self):
        if not self.peaks.any():
            return 0

        return self.peaks[np.argmax(self.__peakHeights)]

    @property
    def __peakHeights(self) -> List[float]:
        return self.peakProperties["peak_heights"]
