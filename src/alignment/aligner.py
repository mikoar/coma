from __future__ import annotations

import itertools
from itertools import dropwhile, takewhile, chain
from typing import List, NamedTuple

from src.alignment.alignment_position import AlignedPair, NotAlignedQueryPosition, \
    NotAlignedReferencePosition, NotAlignedPosition, AlignmentPosition
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentsWithResolvedConflicts
from src.alignment.segments_factory import AlignmentSegmentsFactory
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

    def align(self, reference: OpticalMap, query: OpticalMap, peakPositions: int | List[int],
              isReverse: bool = False) -> AlignmentResultRow:
        if isinstance(peakPositions, int):
            peakPositions = [peakPositions]
        segments = list(itertools.chain.from_iterable(
            [self.getSegments(isReverse, p, query, reference) for p in peakPositions]))

        return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts.create(segments),
                                         query.moleculeId,
                                         reference.moleculeId,
                                         query.length,
                                         reference.length,
                                         isReverse)

    def getSegments(self, isReverse, peakPosition, query, reference):
        referenceStartPosition = peakPosition
        referenceEndPosition = peakPosition + query.length
        alignmentPositions = self.alignmentEngine.align(reference, query, referenceStartPosition,
                                                        referenceEndPosition, isReverse)
        scoredPositions = self.scorer.getScoredPositions(alignmentPositions)
        return self.segmentsFactory.getSegments(scoredPositions)
