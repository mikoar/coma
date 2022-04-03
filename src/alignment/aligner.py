from __future__ import annotations

from itertools import dropwhile, takewhile, chain
from typing import List, NamedTuple

from src.alignment.alignment_position import AlignedPair, nullAlignedPair, NotAlignedQueryPosition, \
    NotAlignedReferencePosition, NotAlignedPosition, AlignmentPosition
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segments import AlignmentSegmentsFactory
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


class AlignerEngine:
    def __init__(self, maxDistance: int):
        self.maxDistance = maxDistance
        self.iteration = 1

    def align(self, reference: OpticalMap, query: OpticalMap, referenceStartPosition: int, referenceEndPosition: int,
              isReverse: bool) -> List[AlignmentPosition]:
        referencePositions = self.__getReferencePositionsWithinRange(reference, referenceStartPosition,
                                                                     referenceEndPosition)
        queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition)
        deduplicatedAlignedPairs = list(AlignedPair.deduplicate(alignedPairs))
        notAlignedPositions = self.__getNotAlignedPositions(queryPositions, referencePositions,
                                                            deduplicatedAlignedPairs, referenceStartPosition)
        return sorted(chain(deduplicatedAlignedPairs, notAlignedPositions))

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
                                  queryPosition.position - referencePositionAdjustedToQuery, self.iteration)
        self.iteration += 1

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


class Aligner:
    def __init__(self, scorer: AlignmentPositionScorer,
                 segmentsFactory: AlignmentSegmentsFactory, alignmentEngine: AlignerEngine) -> None:
        self.scorer = scorer
        self.segmentsFactory = segmentsFactory
        self.alignmentEngine = alignmentEngine

    def align(self, reference: OpticalMap, query: OpticalMap, peakPosition: int,
              isReverse: bool = False) -> AlignmentResultRow:
        referenceStartPosition = peakPosition
        referenceEndPosition = peakPosition + query.length
        alignmentPositions = self.alignmentEngine.align(reference, query, referenceStartPosition, referenceEndPosition,
                                                        isReverse)

        scoredPositions = self.scorer.getScoredPositions(alignmentPositions)
        segments = self.segmentsFactory.getSegments(scoredPositions)

        # TODO: refactor AlignmentResultRow logic
        filteredPositions = sorted(p for s in segments for p in s.positions if isinstance(p, AlignedPair))

        firstPair = filteredPositions[0] if filteredPositions else nullAlignedPair
        lastPair = filteredPositions[-1] if filteredPositions else nullAlignedPair
        queryStart = query.positions[firstPair.query.siteId - 1] if query.positions else 0
        queryEnd = query.positions[lastPair.query.siteId - 1] if query.positions else 0

        return AlignmentResultRow(segments,
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
