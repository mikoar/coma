from bisect import bisect_left
from collections import namedtuple
from dataclasses import dataclass
from itertools import dropwhile, takewhile
from typing import Iterator, List

from .optical_map import OpticalMap, PositionWithSiteId

AlignedPair = namedtuple("AlignedPair", ["reference", "query"])


@dataclass
class AlignmentResult:
    referenceStartPosition: int
    referenceEndPosition: int
    alignedPairs: List[AlignedPair]


class Aligner:
    def __init__(self, maxDistance: int) -> None:
        self.maxDistance = maxDistance

    def align(self, reference: OpticalMap, query: OpticalMap, peakPosition: int, reverse: bool = False) -> AlignmentResult:
        referenceStartPosition = round(peakPosition - query.length/2)
        referenceEndPosition = round(peakPosition + query.length/2)

        referencePositions = list(takewhile(lambda x: x.position <= referenceEndPosition, dropwhile(
            lambda x: x.position < referenceStartPosition - self.maxDistance, reference.getPositionsWithSiteIds())))

        queryPositions = query.getPositionsWithSiteIds(reverse)

        alignedPairs = list(self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition))

        return AlignmentResult(referenceStartPosition, referenceEndPosition, alignedPairs)

    def __getAlignedPairs(self, referencePositions: List[PositionWithSiteId],  queryPositions: Iterator[PositionWithSiteId], shift: int):
        referencePositionsWithoutIndices = list(map(lambda x: x.position, referencePositions))
        for queryPosition in queryPositions:
            referenceIndex = self.__findClosestPositionInReference(referencePositionsWithoutIndices, queryPosition.position + shift)
            if referenceIndex != None:
                referencePosition = referencePositions[referenceIndex]
                yield AlignedPair(referencePosition.siteId, queryPosition.siteId)

    def __findClosestPositionInReference(self, referencePositions: List[int], queryPosition: int):
        index = bisect_left(referencePositions, queryPosition)

        if index == 0:
            return self.__returnIfisWithinMaxDistance(index, referencePositions[index], queryPosition)
        if index == len(referencePositions):
            return self.__returnIfisWithinMaxDistance(-1, referencePositions[-1], queryPosition)
        before = referencePositions[index - 1]
        after = referencePositions[index]
        beforeDistance = queryPosition - before
        afterDistance = after - queryPosition
        if afterDistance < beforeDistance:
            return self.__returnIfisWithinMaxDistance(index, after, queryPosition)
        else:
            return self.__returnIfisWithinMaxDistance(index - 1, before, queryPosition)

    def __returnIfisWithinMaxDistance(self, returned: int, referencePosition: int, queryPosition: int):
        return returned if abs(referencePosition - queryPosition) <= self.maxDistance else None
