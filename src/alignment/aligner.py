from __future__ import annotations

from itertools import dropwhile, takewhile, groupby
from typing import Iterator, List, NamedTuple

from src.alignment.aligned_pair import AlignedPair, nullAlignedPair, NotAlignedQueryPosition
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
        referenceWindowStartPosition = peakPosition
        referenceWindowEndPosition = peakPosition + query.length

        referencePositions = self.__getReferencePositionsWithinRange(reference, referenceWindowStartPosition,
                                                                     referenceWindowEndPosition)

        queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceWindowStartPosition)
        deduplicatedAlignedPairs = list(self.__removeDuplicatedReferenceAlignments(alignedPairs))

        firstPair = deduplicatedAlignedPairs[0] if deduplicatedAlignedPairs else nullAlignedPair
        lastPair = deduplicatedAlignedPairs[-1] if deduplicatedAlignedPairs else nullAlignedPair
        queryStart = query.positions[firstPair.queryPositionIndex - 1]
        queryEnd = query.positions[lastPair.queryPositionIndex - 1]

        notAlignedPositions = self.__getNotAlignedPositions(queryPositions, deduplicatedAlignedPairs)
        return AlignmentResultRow(self.__concatAndSort(deduplicatedAlignedPairs, notAlignedPositions, isReverse),
                                  query.moleculeId,
                                  reference.moleculeId,
                                  *((queryEnd, queryStart) if isReverse else (
                                      queryStart, queryEnd)),
                                  reference.positions[firstPair.referencePositionIndex - 1],
                                  reference.positions[lastPair.referencePositionIndex - 1],
                                  query.length,
                                  reference.length,
                                  isReverse)

    def __getReferencePositionsWithinRange(self, reference: OpticalMap, referenceStartPosition: int,
                                           referenceEndPosition: int):
        return list(takewhile(lambda x: x.position <= referenceEndPosition + self.maxDistance, dropwhile(
            lambda x: x.position < referenceStartPosition - self.maxDistance,
            reference.getPositionsWithSiteIds())))

    def __getAlignedPairs(self, referencePositions: List[PositionWithSiteId],
                          queryPositions: List[PositionWithSiteId], referenceStartPosition: int):
        for referencePosition in referencePositions:
            referencePositionAdjustedToQuery = referencePosition.position - referenceStartPosition
            queryPositionsWithinDistance = takewhile(
                lambda x: x.position <= referencePositionAdjustedToQuery + self.maxDistance, dropwhile(
                    lambda x: x.position < referencePositionAdjustedToQuery - self.maxDistance,
                    queryPositions))
            for queryPosition in queryPositionsWithinDistance:
                yield AlignedPair(referencePosition.siteId, queryPosition.siteId,
                                  queryPosition.position - referencePositionAdjustedToQuery)

    @staticmethod
    def __removeDuplicatedReferenceAlignments(pairs: Iterator[AlignedPair]):
        for _, ambiguousPairs in groupby(pairs, AlignedPair.referenceSelector):
            yield min(ambiguousPairs, key=AlignedPair.distanceSelector)

    @staticmethod
    def __getNotAlignedPositions(queryPositions: List[PositionWithSiteId],
                                 alignedPairs: List[AlignedPair]):
        alignedSiteIds = [p.queryPositionIndex for p in alignedPairs]
        return [NotAlignedQueryPosition(q.siteId) for q in queryPositions if q.siteId not in alignedSiteIds]

    @staticmethod
    def __concatAndSort(l1: List[NotAlignedQueryPosition | AlignedPair],
                        l2: List[NotAlignedQueryPosition | AlignedPair],
                        isReverse: bool):
        return sorted(l1 + l2, key=AlignedPair.querySelector, reverse=isReverse)
