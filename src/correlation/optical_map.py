from __future__ import annotations

from dataclasses import dataclass
from math import ceil
from typing import List

import numpy as np
from scipy.signal import find_peaks, correlate

from src.correlation.bionano_alignment import BionanoAlignment
from src.correlation.peak import Peak
from src.correlation.sequence_generator import SequenceGenerator
from src.diagnostic.validator import Validator


@dataclass(frozen=True)
class PositionWithSiteId:
    siteId: int
    position: int

    def __lt__(self, other: PositionWithSiteId):
        return self.position < other.position


def adjustPeakPositions(peakPositions: np.ndarray, resolution: int, start: int = 0) -> np.ndarray:
    resolutionAdjustment = ceil(resolution / 2) - 1
    return peakPositions * resolution + (resolutionAdjustment + start)


@dataclass(frozen=True)
class OpticalMap:
    moleculeId: int
    length: int
    positions: List[int]

    def trim(self):
        if not self.positions:
            return self
        return OpticalMap(self.moleculeId,
                          self.positions[-1] - self.positions[0] + 1,
                          list(map(lambda p: p - self.positions[0], self.positions)))

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
        if self.length > reference.length:
            return EmptyInitialAlignment(self, reference, sequenceGenerator.resolution, sequenceGenerator.blurRadius)

        sequence = self.getSequence(sequenceGenerator, reverseStrand)
        referenceSequence = reference.getSequence(sequenceGenerator)
        correlation = self.__getCorrelation(referenceSequence, sequence)
        if flatten:
            normalizingFactor = self.__getCorrelation(referenceSequence, np.ones(len(sequence))) + np.sum(sequence)
            correlation = correlation / normalizingFactor

        correlation /= np.max(correlation)

        peakPositions, peakProperties = find_peaks(correlation, height=0.75, width=(None, None), rel_height=0.5,
                                                   distance=((5 * 10 ** 6) / sequenceGenerator.resolution))
        return InitialAlignment(correlation, self, reference, peakPositions, peakProperties, reverseStrand,
                                sequenceGenerator.resolution, sequenceGenerator.blurRadius, 0,
                                len(correlation) * sequenceGenerator.resolution)

    def getSequence(self, sequenceGenerator: SequenceGenerator, reverseStrand=False, start: int = 0, end: int = None):
        sequence = sequenceGenerator.positionsToSequence(self.positions, start, end)
        return sequence[::-1] if reverseStrand else sequence

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='valid', method='fft')


@dataclass()
class CorrelationResult:
    correlation: np.ndarray
    query: OpticalMap
    reference: OpticalMap
    peakPositions: np.ndarray
    peakProperties: dict
    reverseStrand: bool
    resolution: int = 1
    blur: int = 0
    correlationStart: int = 0
    correlationEnd: int = None

    def __post_init__(self):
        self.peakPositions = adjustPeakPositions(self.peakPositions, self.resolution, self.correlationStart)
        self.peakProperties["left_ips"] = adjustPeakPositions(self.peakProperties["left_ips"], self.resolution,
                                                              self.correlationStart)
        self.peakProperties["right_ips"] = adjustPeakPositions(self.peakProperties["right_ips"], self.resolution,
                                                               self.correlationStart)

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
        return self.__getScoreRelativeToNthPeak(peakHeight)

    def getScore(self):
        return self.maxPeak.height - self.correlation.mean() if self.maxPeak else 0

    @property
    def maxPeak(self):
        if not self.peakPositions.any():
            return

        maxIndex = np.argmax(self.__peakHeights)
        return Peak(self.peakPositions[maxIndex], self.__peakHeights[maxIndex], self.__leftBases[maxIndex],
                    self.__rightBases[maxIndex])

    @property
    def peaks(self):
        for position, height, leftBase, rightBase in zip(self.peakPositions, self.__peakHeights, self.__leftBases,
                                                         self.__rightBases):
            yield Peak(position, height, leftBase, rightBase)

    @property
    def __peakHeights(self) -> np.ndarray:
        return self.peakProperties["peak_heights"]

    @property
    def __leftBases(self) -> np.ndarray:
        return self.peakProperties["left_ips"]

    @property
    def __rightBases(self) -> np.ndarray:
        return self.peakProperties["right_ips"]

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

    def __getScoreRelativeToNthPeak(self, peakHeight: float = None):
        heights = self.__peakHeights
        if len(heights) == 0:
            return 0.

        if len(heights) == 1:
            return 1.

        order = min(5, len(heights))

        firstPeakHeight, nthPeakHeight = sorted(heights, reverse=True)[:order:order - 1]
        peakHeightToScore = peakHeight if peakHeight else firstPeakHeight
        return peakHeightToScore - nthPeakHeight


class InitialAlignment(CorrelationResult):
    def refine(self, sequenceGenerator: SequenceGenerator, secondaryMargin: int = 8000,
               peakHeightThreshold: float = 15.):
        querySequence = self.query.getSequence(sequenceGenerator, self.reverseStrand)
        peakPosition = self.maxPeak.position
        resolution = sequenceGenerator.resolution
        referenceStart = peakPosition - secondaryMargin
        referenceEnd = peakPosition + self.query.length + secondaryMargin
        referenceSequence = self.reference.getSequence(sequenceGenerator, False, referenceStart, referenceEnd)
        correlation = self.__getCorrelation(referenceSequence, querySequence)
        peakPositions, peakProperties = find_peaks(correlation, height=peakHeightThreshold,
                                                   width=(None, None), prominence=0.05 * correlation.max())

        correlationLength = len(correlation) * resolution

        return CorrelationResult(correlation, self.query, self.reference, peakPositions, peakProperties,
                                 self.reverseStrand, resolution, sequenceGenerator.blurRadius,
                                 referenceStart,
                                 referenceStart + correlationLength)

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='valid', method='fft')


class EmptyInitialAlignment(InitialAlignment):
    def __init__(self,
                 query: OpticalMap,
                 reference: OpticalMap,
                 resolution: int,
                 blur: int):
        self.query = query
        self.reference = reference
        self.resolution = resolution
        self.blur = blur

    def getScore(self):
        return 0.

    @property
    def peaks(self):
        return []
