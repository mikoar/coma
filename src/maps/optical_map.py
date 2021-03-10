from __future__ import annotations
from typing import List, Literal
from maps.reference_optical_map import ReferenceOpticalMap
from processing.correlate import crossCorrelate


class CorrelationResult(object):
    def __init__(self, correlation, query: OpticalMap, reference: OpticalMap | ReferenceOpticalMap) -> None:
        self.correlation = correlation
        self.query = query
        self.reference = reference


class PositionRange(object):
    def __init__(self, start: int, stop: int) -> None:
        self.start = start
        self.stop = stop


class OpticalMap(ReferenceOpticalMap):
    def __init__(self, moleculeId, strand: Literal["+", "-", "unknown"], sequence: List[int],
                 referenceCoordinates: PositionRange) -> None:
        if strand == "-":
            sequence.reverse()

        super().__init__(sequence)
        self.moleculeId = moleculeId
        self.strand = strand
        self.referenceCoordinates = referenceCoordinates

    def correlate(self, reference: OpticalMap | ReferenceOpticalMap):
        correlation = crossCorrelate(reference.sequence, self.sequence)
        return CorrelationResult(correlation, self, reference)
