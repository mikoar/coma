from __future__ import annotations
from typing import List
import numpy as np

from typing import Iterable
from scipy.signal import correlate
from scipy.signal import find_peaks


class ReferenceOpticalMap:
    def __init__(self, moleculeId: int, sequence:  List[int], positions:  List[int], resolution: int) -> None:
        self.sequence = sequence
        self.positions = positions
        self.moleculeId = moleculeId
        self.resolution = resolution
        self.reverseStrand = False


class OpticalMap:
    def __init__(self, moleculeId: int, sequence: List[int], positions:  List[int], resolution: int) -> None:
        self.sequence = sequence
        self.positions = positions
        self.moleculeId = moleculeId
        self.resolution = resolution
        self.reverseStrand = False  # TODO

    def reverse(self):
        self.sequence.reverse()

    def correlate(self, reference: ReferenceOpticalMap, flatten=True):
        correlation = self.__getCorrelation(reference.sequence, self.sequence)

        if flatten:
            normalizingFactor = self.__getCorrelation(reference.sequence, [1] * len(self.sequence)) + sum(self.sequence)
            correlation = correlation / normalizingFactor

        correlation /= np.max(correlation)

        peaks, peakProperties = find_peaks(correlation, width=(0, len(self.sequence)), height=0.7, distance=(10 ** 5 / self.resolution))
        return CorrelationResult(correlation, self, reference, peaks, peakProperties)

    def __getCorrelation(self, reference: Iterable[int], query: Iterable[int]):
        return correlate(reference, query, mode='same')


class CorrelationResult:
    def __init__(self, correlation, query: OpticalMap, reference: ReferenceOpticalMap, peaks: List[int], peakProperties) -> None:
        self.correlation = correlation
        self.query = query
        self.reference = reference
        self.peaks = peaks
        self.peakProperties = peakProperties
