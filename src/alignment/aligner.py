from __future__ import annotations

from itertools import dropwhile, takewhile, groupby
from typing import Iterator, List, NamedTuple

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
        referenceWindowStartPosition = peakPosition
        referenceWindowEndPosition = peakPosition + query.length

        referencePositions = self.__getReferencePositionsWithinRange(reference, referenceWindowStartPosition,
                                                                     referenceWindowEndPosition)

        queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceWindowStartPosition)
        deduplicatedAlignedPairs = list(self.__removeDuplicatedReferenceAlignments(alignedPairs))

        queryStart = query.positions[deduplicatedAlignedPairs[0].queryPositionIndex - 1]
        queryEnd = query.positions[deduplicatedAlignedPairs[-1].queryPositionIndex - 1]

        return AlignmentResultRow(deduplicatedAlignedPairs,
                                  query.moleculeId,
                                  reference.moleculeId,
                                  *((queryEnd, queryStart) if isReverse else (
                                      queryStart, queryEnd)),
                                  reference.positions[deduplicatedAlignedPairs[0].referencePositionIndex - 1],
                                  reference.positions[deduplicatedAlignedPairs[-1].referencePositionIndex - 1],
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
