from __future__ import annotations
from bisect import bisect_left
from collections import namedtuple
from dataclasses import dataclass
from itertools import dropwhile, takewhile, groupby
from typing import Callable, Iterator, List, NamedTuple, Tuple

from .optical_map import OpticalMap, PositionWithSiteId


class AlignedPair(NamedTuple):
    referencePositionIndex: int
    queryPositionIndex: int
    queryShift: int = 0

    @property
    def distance(self):
        return abs(self.queryShift)

    def getFalseNegativesCount(self, previousPair: AlignedPair | None):
        return 0 if not previousPair else self.referencePositionIndex - previousPair.referencePositionIndex - 1

    def getFalsePositivesCount(self, previousPair: AlignedPair | None):
        return 0 if not previousPair else abs(self.queryPositionIndex - previousPair.queryPositionIndex) - 1

    def __repr__(self) -> str:
        return f"({self.referencePositionIndex}, {self.queryPositionIndex})"

    def __eq__(self, other) -> bool:
        return self.referencePositionIndex == other[0] and self.queryPositionIndex == other[1]


@dataclass
class AlignmentResult:
    referenceStartPosition: int
    referenceEndPosition: int
    alignedPairs: List[AlignedPair]


class _ReferenceIndexWithDistance(NamedTuple):
    index: int
    distance: int
    queryShift: int

    @staticmethod
    def withQueryAfterReference(index: int, distance: int):
        return _ReferenceIndexWithDistance(index, distance, distance)

    @staticmethod
    def withQueryBeforeReference(index: int, distance: int):
        return _ReferenceIndexWithDistance(index, distance, -distance)


class Aligner:

    def __init__(self, maxDistance: int) -> None:
        self.maxDistance = maxDistance

    def align(self, reference: OpticalMap, query: OpticalMap, peakPosition: int, isReverse: bool = False) -> AlignmentResult:
        referenceStartPosition = round(peakPosition - query.length/2)
        referenceEndPosition = round(peakPosition + query.length/2)

        referencePositions = list(takewhile(lambda x: x.position <= referenceEndPosition + self.maxDistance, dropwhile(
            lambda x: x.position < referenceStartPosition - self.maxDistance, reference.getPositionsWithSiteIds())))

        queryPositions = query.getPositionsWithSiteIds(isReverse)
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition)
        deduplicatedAlignedPairs = self.__removeDuplicatedPairsWithNonMinimalDistance(alignedPairs)

        return AlignmentResult(referenceStartPosition, referenceEndPosition, list(deduplicatedAlignedPairs))

    def __getAlignedPairs(self, referencePositions: List[PositionWithSiteId],
                          queryPositions: Iterator[PositionWithSiteId], referenceStartPosition: int):
        referencePositionsWithoutIndices = list(map(lambda x: x.position, referencePositions))
        for queryPosition in queryPositions:
            indexWithDistance = self.__findClosestPositionInReference(referencePositionsWithoutIndices,
                                                                      queryPosition.position + referenceStartPosition)
            if indexWithDistance is not None:
                referencePosition = referencePositions[indexWithDistance.index]
                yield AlignedPair(referencePosition.siteId, queryPosition.siteId, indexWithDistance.queryShift)

    def __findClosestPositionInReference(self, referencePositions: List[int], queryPosition: int):
        index = bisect_left(referencePositions, queryPosition)

        if index == 0:
            distanceToNextReferencePosition = self.__getDistanceToNextReferencePosition(queryPosition, referencePositions, index)
            return self.__returnIfisWithinMaxDistance(_ReferenceIndexWithDistance.withQueryBeforeReference(index, distanceToNextReferencePosition))
        if index == len(referencePositions):
            distanceToPreviousReferencePosition = self.__getDistanceToPreviousReferencePosition(queryPosition, referencePositions, index)
            return self.__returnIfisWithinMaxDistance(_ReferenceIndexWithDistance.withQueryAfterReference(-1, distanceToPreviousReferencePosition))

        distanceToNextReferencePosition = self.__getDistanceToNextReferencePosition(queryPosition, referencePositions, index)
        distanceToPreviousReferencePosition = self.__getDistanceToPreviousReferencePosition(queryPosition, referencePositions, index)
        if distanceToNextReferencePosition < distanceToPreviousReferencePosition:
            return self.__returnIfisWithinMaxDistance(_ReferenceIndexWithDistance.withQueryBeforeReference(index, distanceToNextReferencePosition))
        else:
            return self.__returnIfisWithinMaxDistance(_ReferenceIndexWithDistance.withQueryAfterReference(index - 1, distanceToPreviousReferencePosition))

    def __getDistanceToNextReferencePosition(self, queryPosition, referencePositions, index):
        nextRef = referencePositions[index]
        return nextRef - queryPosition

    def __getDistanceToPreviousReferencePosition(self, queryPosition, referencePositions, index):
        previousRefPosition = referencePositions[index - 1]
        return queryPosition - previousRefPosition

    def __returnIfisWithinMaxDistance(self, pair: _ReferenceIndexWithDistance):
        return pair if pair.distance <= self.maxDistance else None

    def __removeDuplicatedPairsWithNonMinimalDistance(self, pairs: Iterator[AlignedPair]):
        referenceSelector: Callable[[AlignedPair], int] = lambda pair: pair.referencePositionIndex
        distanceSelector: Callable[[AlignedPair], int] = lambda pair: pair.distance
        for _, ambiguousPairs in groupby(pairs, referenceSelector):
            yield min(ambiguousPairs, key=distanceSelector)
