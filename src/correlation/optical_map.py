from __future__ import annotations

import warnings
from dataclasses import dataclass
from math import ceil
from typing import List

import numpy as np
from scipy.signal import find_peaks, correlate

from src.correlation.peak import Peak
from src.correlation.sequence_generator import SequenceGenerator

warnings.simplefilter("ignore")


@dataclass(frozen=True)
class PositionWithSiteId:
    siteId: int
    position: int

    def __lt__(self, other: PositionWithSiteId):
        return self.position < other.position


def toRelativeGenomicPositions(correlationCoordinates: np.ndarray, resolution: int, start: int = 0) -> np.ndarray:
    resolutionAdjustment = ceil(resolution / 2) - 1
    return correlationCoordinates * resolution + (resolutionAdjustment + start)


@dataclass(frozen=True)
class OpticalMap:
    moleculeId: int
    length: int
    positions: List[int]
    shift: int = 0

    def trim(self):
        if not self.positions:
            return self
        return OpticalMap(self.moleculeId,
                          self.positions[-1] - self.positions[0] + 1,
                          list(map(lambda p: p - self.positions[0], self.positions)))

    def getPositionsWithSiteIds(self, reverse: bool = False):
        if reverse:
            i = len(self.positions) + self.shift
            moleculeEndPosition = self.length - 1
            for position in self.positions[::-1]:
                yield PositionWithSiteId(i, moleculeEndPosition - position)
                i -= 1
        else:
            i = 1 + self.shift
            for position in self.positions:
                yield PositionWithSiteId(i, position)
                i += 1

    def getInitialAlignment(self, reference: OpticalMap, sequenceGenerator: SequenceGenerator, minPeakDistance: int,
                            peaksCount: int, reverseStrand=False):
        if self.length > reference.length:
            return EmptyInitialAlignment(self, reference, sequenceGenerator.resolution, sequenceGenerator.blurRadius)

        sequence = self.getSequence(sequenceGenerator, reverseStrand)
        referenceSequence = reference.getSequence(sequenceGenerator)
        correlation = self.__getCorrelation(referenceSequence, sequence)

        normalizingFactor = (self.__getCorrelation(referenceSequence, np.ones(len(sequence))) + np.sum(sequence)) / 2
        correlation = correlation / normalizingFactor

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            peakPositions, peakProperties = find_peaks(
                correlation,
                height=0.75 * np.max(correlation),
                width=(None, None),
                rel_height=0.5,
                distance=(minPeakDistance / sequenceGenerator.resolution))

        return InitialAlignment.create(correlation, self, reference, peakPositions, peakProperties, peaksCount,
                                       reverseStrand,
                                       sequenceGenerator.resolution, sequenceGenerator.blurRadius, 0,
                                       len(correlation) * sequenceGenerator.resolution)

    def getSequence(self, sequenceGenerator: SequenceGenerator, reverseStrand=False, start: int = 0, end: int = None):
        sequence = sequenceGenerator.positionsToSequence(self.positions, start, end)
        return sequence[::-1] if reverseStrand else sequence

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='valid', method='fft')


class CorrelationResult:
    @staticmethod
    def create(correlation: np.ndarray,
               query: OpticalMap,
               reference: OpticalMap,
               peakPositions: np.ndarray,
               peakProperties: dict,
               peaksCount: int,
               reverseStrand: bool,
               resolution: int = 1,
               blur: int = 0,
               correlationStart: int = 0,
               correlationEnd: int = None,
               peakHeightThreshold: float = None):
        return CorrelationResult(
            correlation,
            query,
            reference,
            CorrelationResult.createPeaks(peakPositions, peakProperties, resolution, correlationStart, 0, peaksCount),
            reverseStrand,
            peakHeightThreshold,
            resolution,
            blur,
            correlationStart,
            correlationEnd or len(correlation) - 1)

    def __init__(self,
                 correlation: np.ndarray,
                 query: OpticalMap,
                 reference: OpticalMap,
                 peaks: List[Peak],
                 reverseStrand: bool,
                 noiseLevel: float,
                 resolution: int = 1,
                 blur: int = 0,
                 correlationStart: int = 0,
                 correlationEnd: int = None):
        self.correlation = correlation
        self.query = query
        self.reference = reference
        self.peaks = peaks
        self.reverseStrand = reverseStrand
        self.peakBaseLevel = noiseLevel
        self.resolution = resolution
        self.blur = blur
        self.correlationStart = correlationStart
        self.correlationEnd = correlationEnd

    @staticmethod
    def createPeaks(peakPositions: np.ndarray, peakProperties: dict, resolution: int, correlationStart: int,
                    noiseLevel: float, peaksCount: int):
        if peaksCount < peakPositions.size:
            bestPeaksIndices = np.argpartition(-peakProperties["peak_heights"], peaksCount)[:peaksCount]
        else:
            bestPeaksIndices = np.arange(peakPositions.size)
        return [Peak(position, height, leftBase, rightBase, height - noiseLevel)
                for position, height, leftBase, rightBase
                in zip(toRelativeGenomicPositions(peakPositions[bestPeaksIndices], resolution, correlationStart),
                       peakProperties["peak_heights"][bestPeaksIndices],
                       toRelativeGenomicPositions(peakProperties["left_ips"][bestPeaksIndices], resolution,
                                                  correlationStart),
                       toRelativeGenomicPositions(peakProperties["right_ips"][bestPeaksIndices], resolution,
                                                  correlationStart))]

    @staticmethod
    def rootMeanSquare(array: np.ndarray) -> float:
        return np.sqrt(np.mean(array[array != 0] ** 2))

    def getScore(self):
        return self.maxPeak.score if self.maxPeak else 0

    @property
    def maxPeak(self):
        return max(self.peaks, key=lambda p: p.height, default=None)


class InitialAlignment(CorrelationResult):
    @staticmethod
    def create(correlation: np.ndarray,
               query: OpticalMap,
               reference: OpticalMap,
               peakPositions: np.ndarray,
               peakProperties: dict,
               peaksCount: int,
               reverseStrand: bool,
               resolution: int = 1,
               blur: int = 0,
               correlationStart: int = 0,
               correlationEnd: int = None,
               peakHeightThreshold: float = None):
        noiseLevel = InitialAlignment.rootMeanSquare(correlation)
        return InitialAlignment(
            correlation,
            query,
            reference,
            InitialAlignment.createPeaks(peakPositions, peakProperties, resolution, correlationStart, noiseLevel,
                                         peaksCount),
            reverseStrand,
            noiseLevel,
            resolution,
            blur,
            correlationStart,
            correlationEnd or len(correlation) - 1)

    def refine(self, peakPosition: int, sequenceGenerator: SequenceGenerator, secondaryMargin: int = 8000,
               peakHeightThreshold: float = 15.):
        querySequence = self.query.getSequence(sequenceGenerator, self.reverseStrand)
        resolution = sequenceGenerator.resolution
        referenceStart = peakPosition - secondaryMargin
        referenceEnd = peakPosition + self.query.length + secondaryMargin
        referenceSequence = self.reference.getSequence(sequenceGenerator, False, referenceStart, referenceEnd)
        correlation = self.__getCorrelation(referenceSequence, querySequence)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            peakPositions, peakProperties = find_peaks(
                correlation,
                height=peakHeightThreshold,
                width=(None, None),
                prominence=0.05 * correlation.max(initial=0))

        correlationLength = len(correlation) * resolution

        return CorrelationResult.create(correlation, self.query, self.reference, peakPositions, peakProperties,
                                        10, self.reverseStrand, resolution, sequenceGenerator.blurRadius,
                                        referenceStart, referenceStart + correlationLength, peakHeightThreshold)

    @staticmethod
    def __getCorrelation(reference: np.ndarray, query: np.ndarray) -> np.ndarray:
        return correlate(reference, query, mode='valid', method='fft')


class EmptyInitialAlignment(InitialAlignment):
    def __init__(self,
                 query: OpticalMap,
                 reference: OpticalMap,
                 resolution: int,
                 blur: int):
        super().__init__(np.array([]),
                         query,
                         reference,
                         [],
                         False,
                         0.)
        self.query = query
        self.reference = reference
        self.resolution = resolution
        self.blur = blur

    def getScore(self):
        return 0.
