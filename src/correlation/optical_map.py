from __future__ import annotations

from dataclasses import dataclass
from typing import List, NamedTuple

import numpy as np
from scipy.signal import correlate, find_peaks

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
        sequence = self.__getSequence(sequenceGenerator, reverseStrand)
        referenceSequence = reference.__getSequence(sequenceGenerator)
        correlation = self.__getCorrelation(referenceSequence, sequence)
        if flatten:
            normalizingFactor = self.__getCorrelation(referenceSequence, np.ones(len(sequence))) + np.sum(sequence)
            correlation = correlation / normalizingFactor

        correlation /= np.max(correlation)

        peakPositions, peakProperties = find_peaks(
            correlation,
            height=0.01,
            distance=((5 * 10 ** 6) / sequenceGenerator.resolution))

        return InitialAlignment(correlation, self, reference, sequenceGenerator.resolution,
                                sequenceGenerator.blurRadius, peakPositions, peakProperties)

    def __getSequence(self, sequenceGenerator: SequenceGenerator, reverseStrand=False):
        sequence = sequenceGenerator.positionsToSequence(self.positions)
        return sequence[::-1] if reverseStrand else sequence

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='same', method='fft')


@dataclass()
class InitialAlignment:
    correlation: np.ndarray
    query: OpticalMap
    reference: OpticalMap
    resolution: int
    blur: int
    peakPositions: np.ndarray
    peakProperties: dict

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
