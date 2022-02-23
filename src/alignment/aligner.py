from __future__ import annotations

from itertools import dropwhile, takewhile, chain
from typing import List, NamedTuple

from src.alignment.alignment_position import AlignedPair, nullAlignedPair, NotAlignedQueryPosition, \
    NotAlignedReferencePosition, NotAlignedPosition
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
        referenceStartPosition = peakPosition
        referenceEndPosition = peakPosition + query.length

        referencePositions = self.__getReferencePositionsWithinRange(reference, referenceStartPosition,
                                                                     referenceEndPosition)

        queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition)
        deduplicatedAlignedPairs = list(AlignedPair.deduplicate(alignedPairs))

        firstPair = deduplicatedAlignedPairs[0] if deduplicatedAlignedPairs else nullAlignedPair
        lastPair = deduplicatedAlignedPairs[-1] if deduplicatedAlignedPairs else nullAlignedPair
        queryStart = query.positions[firstPair.query.siteId - 1] if query.positions else 0
        queryEnd = query.positions[lastPair.query.siteId - 1] if query.positions else 0

        notAlignedPositions = self.__getNotAlignedPositions(queryPositions, referencePositions,
                                                            deduplicatedAlignedPairs, referenceStartPosition)

        return AlignmentResultRow(sorted(chain(deduplicatedAlignedPairs, notAlignedPositions)),
                                  query.moleculeId,
                                  reference.moleculeId,
                                  *((queryEnd, queryStart) if isReverse else (
                                      queryStart, queryEnd)),
                                  reference.positions[
                                      firstPair.reference.siteId - 1] if reference.positions else 0,
                                  reference.positions[
                                      lastPair.reference.siteId - 1] if reference.positions else 0,
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
                yield AlignedPair(referencePosition, queryPosition,
                                  queryPosition.position - referencePositionAdjustedToQuery)

    @staticmethod
    def __getNotAlignedPositions(queryPositions: List[PositionWithSiteId],
                                 referencePositions: List[PositionWithSiteId],
                                 alignedPairs: List[AlignedPair],
                                 referenceStartPosition: int):
        alignedReferenceSiteIds = [p.reference.siteId for p in alignedPairs]
        alignedQuerySiteIds = [p.query.siteId for p in alignedPairs]
        notAlignedReferencePositions: List[NotAlignedPosition] = \
            [NotAlignedReferencePosition(r) for r in referencePositions if
             r.siteId not in alignedReferenceSiteIds]
        notAlignedQueryPositions = [NotAlignedQueryPosition(q, referenceStartPosition) for q in queryPositions if
                                    q.siteId not in alignedQuerySiteIds]
        return notAlignedReferencePositions + notAlignedQueryPositions
