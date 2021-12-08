from __future__ import annotations

from dataclasses import dataclass
from math import floor
from typing import List, NamedTuple

import numpy as np
from scipy.signal import find_peaks, correlate

from src.correlation.bionano_alignment import BionanoAlignment
from src.correlation.peak import Peak
from src.correlation.sequence_generator import SequenceGenerator
from src.correlation.validator import Validator


class PositionWithSiteId(NamedTuple):
    siteId: int
    position: int


@dataclass(frozen=True)
class OpticalMap:
    moleculeId: int
    length: int
    positions: List[int]

    def getPositionsWithSiteIds(self, reverse: bool = False):
        if reverse:
            i = len(self.positions)
            moleculeEndPosition = self.length - 1
            for position in self.positions[::-1]:
                yield PositionWithSiteId(i, moleculeEndPosition - position)
                i -= 1
        else:
            i = 1
            for position in self.positions:
                yield PositionWithSiteId(i, position)
                i += 1

    def getInitialAlignment(self, reference: OpticalMap, sequenceGenerator: SequenceGenerator, reverseStrand=False,
                            flatten=True):
        sequence = self.getSequence(sequenceGenerator, reverseStrand)
        referenceSequence = reference.getSequence(sequenceGenerator)
        correlation = self.__getCorrelation(referenceSequence, sequence)
        if flatten:
            normalizingFactor = self.__getCorrelation(referenceSequence, np.ones(len(sequence))) + np.sum(sequence)
            correlation = correlation / normalizingFactor

        correlation /= np.max(correlation)

        peakPositions, peakProperties = find_peaks(correlation, height=0.01,
                                                   distance=((5 * 10 ** 6) / sequenceGenerator.resolution))

        return InitialAlignment(correlation, self, reference, peakPositions, peakProperties, reverseStrand,
                                sequenceGenerator.resolution, sequenceGenerator.blurRadius, 0,
                                len(correlation) * sequenceGenerator.resolution)

    def getSequence(self, sequenceGenerator: SequenceGenerator, reverseStrand=False, start: int = 0, end: int = None):
        sequence = sequenceGenerator.positionsToSequence(self.positions, start, end)
        return sequence[::-1] if reverseStrand else sequence

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='same', method='fft')


@dataclass()
class CorrelationResult:
    correlation: np.ndarray
    query: OpticalMap
    reference: OpticalMap
    peakPositions: np.ndarray
    peakProperties: dict
    reverseStrand: bool
    resolution: int = 1
    blur: int = 1
    correlationStart: int = 0
    correlationEnd: int = None

    def __post_init__(self):
        if self.correlationEnd is None:
            self.correlationEnd = len(self.correlation) - 1

    def getRelativeScore(self, reference: BionanoAlignment, validator: Validator):
        maxPeak = self.maxPeak
        isMaxValid = validator.validate(maxPeak, reference)
        if maxPeak and isMaxValid:
            peakHeight = maxPeak.height
        else:
            peakHeight = self.__getMaxValidPeakHeight(reference,
                                                      validator) or self.__getMaxCorrelationValueInAlignmentRange(
                reference)
        return self.__score(peakHeight)

    @property
    def maxPeak(self):
        if not self.peakPositions.any():
            return

        maxIndex = np.argmax(self.__peakHeights)
        return Peak(self.peakPositions[maxIndex], self.__peakHeights[maxIndex], self.resolution)

    @property
    def peaks(self):
        for position, height in zip(self.peakPositions, self.__peakHeights):
            yield Peak(position, height, self.resolution)

    @property
    def __peakHeights(self) -> List[float]:
        return self.peakProperties["peak_heights"]

    def __getMaxValidPeakHeight(self, reference: BionanoAlignment, validator: Validator):
        validPeaks = [peak for peak in self.peaks if validator.validate(peak, reference)]
        maxValidPeakHeight = max(validPeaks, key=lambda peak: peak.height).height if validPeaks else None
        return maxValidPeakHeight

    def __getMaxCorrelationValueInAlignmentRange(self, reference: BionanoAlignment) -> float:
        expectedQueryStartPosition = int(reference.expectedQueryMoleculeStart / self.resolution)
        expectedQueryEndPosition = int(reference.expectedQueryMoleculeStart / self.resolution)
        expectedQueryRange: np.ndarray = self.correlation[
                                         expectedQueryStartPosition: expectedQueryEndPosition]
        return np.max(expectedQueryRange) if expectedQueryRange.any() else 0.

    def __score(self, peakHeight: float | None):
        order = 5
        heights = self.__peakHeights

        if len(heights) < order:
            return 1.

        firstPeakHeight, nthPeakHeight = sorted(heights, reverse=True)[:order:order - 1]
        peakHeightToScore = peakHeight if peakHeight else firstPeakHeight
        return peakHeightToScore - nthPeakHeight


class InitialAlignment(CorrelationResult):
    def refine(self, sequenceGenerator: SequenceGenerator, maxAdjustment: int):
        querySequence = self.query.getSequence(sequenceGenerator, self.reverseStrand)
        peakPosition = self.maxPeak.positionInReference
        referenceStart = peakPosition - floor(self.query.length / 2) - maxAdjustment
        referenceEnd = peakPosition + floor(self.query.length / 2) + maxAdjustment
        referenceSequence = self.reference.getSequence(sequenceGenerator, self.reverseStrand, referenceStart,
                                                       referenceEnd)
        correlation = self.__getCorrelation(referenceSequence, querySequence)
        peakPositions, peakProperties = find_peaks(correlation, height=10 * sequenceGenerator.blurRadius)

        querySequenceLength = len(querySequence)
        adjustedPeakPositions = self.__adjustPeakPositionsToFullReference(peakPositions, referenceStart,
                                                                          querySequenceLength,
                                                                          sequenceGenerator.resolution)
        correlationLength = len(correlation) * sequenceGenerator.resolution

        return CorrelationResult(correlation, self.query, self.reference, adjustedPeakPositions, peakProperties,
                                 self.reverseStrand, sequenceGenerator.resolution, sequenceGenerator.blurRadius,
                                 peakPosition - floor(correlationLength / 2),
                                 peakPosition + floor(correlationLength / 2))

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='valid', method='fft')

    @staticmethod
    def __adjustPeakPositionsToFullReference(peakPositions: np.ndarray, referenceStart: int, queryLength: int,
                                             resolution: int):
        return peakPositions + referenceStart / resolution + floor(queryLength / 2)
