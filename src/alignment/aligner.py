from __future__ import annotations

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

        referencePositions = self.__getReferencePositionsWithinRange(reference, referenceStartPosition,
                                                                     referenceEndPosition)

        queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition)
        deduplicatedAlignedPairs = self.__removeDuplicatedReferenceAlignments(alignedPairs)

        return AlignmentResultRow(list(deduplicatedAlignedPairs),
                                  query.moleculeId,
                                  reference.moleculeId,
                                  *((query.positions[-1], query.positions[0]) if isReverse else (
                                      query.positions[0], query.positions[-1])),
                                  referenceStartPosition,
                                  referenceEndPosition,
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
        referenceSelector: Callable[[AlignedPair], int] = lambda pair: pair.referencePositionIndex
        distanceSelector: Callable[[AlignedPair], int] = lambda pair: pair.distance
        for _, ambiguousPairs in groupby(pairs, referenceSelector):
            yield min(ambiguousPairs, key=distanceSelector)
