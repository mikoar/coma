from __future__ import annotations
from typing import List
import numpy as np
from processing.correlate import crossCorrelate


class ReferenceOpticalMap:
    def __init__(self, sequence:  List[int], positions:  List[int]) -> None:
        self.sequence = sequence
        self.positions = positions


class OpticalMap(ReferenceOpticalMap):
    def __init__(self, moleculeId: int, sequence: List[int], positions:  List[int]) -> None:
        super().__init__(sequence, positions)
        self.moleculeId = moleculeId

    def reverse(self):
        self.sequence.reverse()

    def correlate(self, reference: ReferenceOpticalMap, normalize=True):
        correlation = crossCorrelate(reference.sequence, self.sequence)

        if normalize:
            normalizingFactor = crossCorrelate(reference.sequence, [1] * len(self.sequence))
            correlation = np.divide(correlation, normalizingFactor)

        return CorrelationResult(correlation, self, reference)


class CorrelationResult:
    def __init__(self, correlation, query: OpticalMap, reference: ReferenceOpticalMap) -> None:
        self.correlation = correlation
        self.query = query
        self.reference = reference
