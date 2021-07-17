from __future__ import annotations
from scipy.signal import find_peaks
from typing import List
import numpy as np
from scipy.signal import correlate


class OpticalMap:
    def __init__(self, moleculeId: int, sequence: np.ndarray, positions:  List[int], resolution: int) -> None:
        self.sequence = sequence
        self.positions = positions
        self.moleculeId = moleculeId
        self.resolution = resolution

    def correlate(self, reference: np.ndarray, reverseStrand=False, flatten=True):
        if reverseStrand:
            self.sequence = self.sequence[::-1]

        correlation = self.__getCorrelation(reference, self.sequence)

        if flatten:
            normalizingFactor = self.__getCorrelation(reference, np.ones(len(self.sequence))) + np.sum(self.sequence)
            correlation = correlation / normalizingFactor

        correlation /= np.max(correlation)

        return CorrelationResult(correlation, self, reference)

    def __getCorrelation(self, reference: np.ndarray, query: np.ndarray):
        return correlate(reference, query, mode='same', method='fft')


class CorrelationResult:
    def __init__(self, correlation: List[float], query: OpticalMap, reference: np.ndarray) -> None:
        self.correlation = correlation
        self.query = query
        self.reference = reference


class Peaks:
    def __init__(self, correlationResult: CorrelationResult) -> None:
        self.correlationResult = correlationResult
        self.peaks, self.peakProperties = find_peaks(correlationResult.correlation,
                                                     height=0.05,
                                                     prominence=0.2,
                                                     distance=((5 * 10 ** 6) / correlationResult.query.resolution))

    @property
    def score(self):
        heights = self.__peakHeights

        if len(heights) == 1:
            return 1.

        if len(heights) < 2:
            return 0.

        twoHighest = sorted(heights, reverse=True)[:2]
        return twoHighest[0] - twoHighest[1]

    @property
    def reverseScore(self):
        score = self.score
        return 1 / score if score else 1

    @property
    def max(self):
        if not self.peaks.any():
            return 0

        return self.peaks[np.argmax(self.__peakHeights)]

    @property
    def __peakHeights(self) -> List[float]:
        return self.peakProperties["peak_heights"]
