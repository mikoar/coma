from __future__ import annotations

from bisect import bisect_left
from itertools import dropwhile, takewhile, groupby
from typing import Callable, Iterator, List, NamedTuple

from src.alignment.aligned_pair import AlignedPair
from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.optical_map import OpticalMap, PositionWithSiteId


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

    def align(self, reference: OpticalMap, query: OpticalMap, peakPosition: int,
              isReverse: bool = False) -> AlignmentResultRow:
        referenceStartPosition = round(peakPosition - query.length / 2)
        referenceEndPosition = round(peakPosition + query.length / 2)

        referencePositions = list(takewhile(lambda x: x.position <= referenceEndPosition + self.maxDistance, dropwhile(
            lambda x: x.position < referenceStartPosition - self.maxDistance, reference.getPositionsWithSiteIds())))

        queryPositions = query.getPositionsWithSiteIds(isReverse)
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition)
        deduplicatedAlignedPairs = self.__removeDuplicatedPairsWithNonMinimalDistance(alignedPairs)

        return AlignmentResultRow(list(deduplicatedAlignedPairs),
                                  query.moleculeId,
                                  reference.moleculeId,
                                  1,
                                  query.length,
                                  referenceStartPosition,
                                  referenceEndPosition,
                                  query.length,
                                  reference.length,
                                  isReverse)

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
            distanceToNextReferencePosition = self.__getDistanceToNextReferencePosition(queryPosition,
                                                                                        referencePositions, index)
            return self.__returnIfIsWithinMaxDistance(
                _ReferenceIndexWithDistance.withQueryBeforeReference(index, distanceToNextReferencePosition))
        if index == len(referencePositions):
            distanceToPreviousReferencePosition = self.__getDistanceToPreviousReferencePosition(queryPosition,
                                                                                                referencePositions,
                                                                                                index)
            return self.__returnIfIsWithinMaxDistance(
                _ReferenceIndexWithDistance.withQueryAfterReference(-1, distanceToPreviousReferencePosition))

        distanceToNextReferencePosition = self.__getDistanceToNextReferencePosition(queryPosition, referencePositions,
                                                                                    index)
        distanceToPreviousReferencePosition = self.__getDistanceToPreviousReferencePosition(queryPosition,
                                                                                            referencePositions, index)
        if distanceToNextReferencePosition < distanceToPreviousReferencePosition:
            return self.__returnIfIsWithinMaxDistance(
                _ReferenceIndexWithDistance.withQueryBeforeReference(index, distanceToNextReferencePosition))
        else:
            return self.__returnIfIsWithinMaxDistance(
                _ReferenceIndexWithDistance.withQueryAfterReference(index - 1, distanceToPreviousReferencePosition))

    def __getDistanceToNextReferencePosition(self, queryPosition, referencePositions, index):
        nextRef = referencePositions[index]
        return nextRef - queryPosition

    def __getDistanceToPreviousReferencePosition(self, queryPosition, referencePositions, index):
        previousRefPosition = referencePositions[index - 1]
        return queryPosition - previousRefPosition

    def __returnIfIsWithinMaxDistance(self, pair: _ReferenceIndexWithDistance):
        return pair if pair.distance <= self.maxDistance else None

    def __removeDuplicatedPairsWithNonMinimalDistance(self, pairs: Iterator[AlignedPair]):
        referenceSelector: Callable[[AlignedPair], int] = lambda pair: pair.referencePositionIndex
        distanceSelector: Callable[[AlignedPair], int] = lambda pair: pair.distance
        for _, ambiguousPairs in groupby(pairs, referenceSelector):
            yield min(ambiguousPairs, key=distanceSelector)
