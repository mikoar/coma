from __future__ import annotations
from typing import List
import numpy as np

from typing import Iterable
from scipy.signal import correlate

from quality import Quality


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

        return CorrelationResult(correlation, self, reference)

    def __getCorrelation(self, reference: Iterable[int], query: Iterable[int]):
        return correlate(reference, query, mode='same')


class CorrelationResult:
    def __init__(self, correlation: List[float], query: OpticalMap, reference: ReferenceOpticalMap) -> None:
        self.correlation = correlation
        self.query = query
        self.reference = reference

    @property
    def quality(self):
        return Quality(self)
