from bisect import bisect_left
from collections import namedtuple
from dataclasses import dataclass
from itertools import dropwhile, takewhile, groupby
from typing import Callable, Iterator, List, NamedTuple, Tuple

from .optical_map import OpticalMap, PositionWithSiteId


class AlignedPair(NamedTuple):
    reference: int
    query: int
    distance: int

    def __repr__(self) -> str:
        return f"({self.reference}, {self.query})"

    def __eq__(self, other) -> bool:
        return self.reference == other[0] and self.query == other[1]


@dataclass
class AlignmentResult:
    referenceStartPosition: int
    referenceEndPosition: int
    alignedPairs: List[AlignedPair]


class _ReferenceIndexWithDistanceToQuery(NamedTuple):
    index: int
    distance: int


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
                yield AlignedPair(referencePosition.siteId, queryPosition.siteId, indexWithDistance.distance)

    def __findClosestPositionInReference(self, referencePositions: List[int], queryPosition: int):
        index = bisect_left(referencePositions, queryPosition)

        if index == 0:
            afterDistance = self.__getDistanceToNextPosition(queryPosition, referencePositions, index)
            return self.__returnIfisWithinMaxDistance(index, afterDistance)
        if index == len(referencePositions):
            beforeDistance = self.__getDistanceToPreviousPosition(queryPosition, referencePositions, index)
            return self.__returnIfisWithinMaxDistance(-1, beforeDistance)

        afterDistance = self.__getDistanceToNextPosition(queryPosition, referencePositions, index)
        beforeDistance = self.__getDistanceToPreviousPosition(queryPosition, referencePositions, index)
        if afterDistance < beforeDistance:
            return self.__returnIfisWithinMaxDistance(index, afterDistance)
        else:
            return self.__returnIfisWithinMaxDistance(index - 1, beforeDistance)

    def __getDistanceToNextPosition(self, queryPosition, referencePositions, index):
        after = referencePositions[index]
        afterDistance = after - queryPosition
        return afterDistance

    def __getDistanceToPreviousPosition(self, queryPosition, referencePositions, index):
        before = referencePositions[index - 1]
        beforeDistance = queryPosition - before
        return beforeDistance

    def __returnIfisWithinMaxDistance(self, index: int, distance: int):
        return _ReferenceIndexWithDistanceToQuery(index, distance) if distance <= self.maxDistance else None

    def __removeDuplicatedPairsWithNonMinimalDistance(self, pairs: Iterator[AlignedPair]):
        referenceSelector: Callable[[AlignedPair], int] = lambda pair: pair.reference
        distanceSelector: Callable[[AlignedPair], int] = lambda pair: pair.distance
        for _, ambiguousPairs in groupby(pairs, referenceSelector):
            yield min(ambiguousPairs, key=distanceSelector)
