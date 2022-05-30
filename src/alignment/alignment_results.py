from __future__ import annotations

import itertools
from dataclasses import dataclass
from enum import Enum
from typing import List

from src.alignment.alignment_position import AlignedPair, NotAlignedPosition
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentsWithResolvedConflicts
from src.alignment.segments import AlignmentSegment


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"


@dataclass
class AlignmentResults:
    referenceFilePath: str
    queryFilePath: str
    rows: List[AlignmentResultRow]


class AlignmentResultRow:
    @staticmethod
    def create(segmentsWithoutConflicts: AlignmentSegmentsWithResolvedConflicts,
               queryId: int,
               referenceId: int,
               queryLength: int,
               referenceLength: int,
               reverseStrand: bool):
        segments = segmentsWithoutConflicts.segments
        alignedPairs = sorted(p for s in segments for p in s.positions if isinstance(p, AlignedPair))
        firstPair = alignedPairs[0] if alignedPairs else AlignedPair.null
        lastPair = alignedPairs[-1] if alignedPairs else AlignedPair.null
        queryStartPosition = (firstPair if not reverseStrand else lastPair).query.position
        queryEndPosition = (lastPair if not reverseStrand else firstPair).query.position
        referenceStartPosition = firstPair.reference.position
        referenceEndPosition = lastPair.reference.position
        confidence = sum(s.segmentScore for s in segments)
        return AlignmentResultRow(segments, queryId, referenceId, queryLength, referenceLength, queryStartPosition,
                                  queryEndPosition, referenceStartPosition, referenceEndPosition, reverseStrand,
                                  confidence)

    def __init__(self,
                 segments: List[AlignmentSegment],
                 queryId: int = 1,
                 referenceId: int = 1,
                 queryLength: int = 1,
                 referenceLength: int = 1,
                 queryStartPosition: int = 0,
                 queryEndPosition: int = 0,
                 referenceStartPosition: int = 0,
                 referenceEndPosition: int = 0,
                 reverseStrand: bool = False,
                 confidence: float = 0.):

        self.segments = segments
        self.queryId = queryId
        self.referenceId = referenceId
        self.queryLength = queryLength
        self.referenceLength = referenceLength
        self.reverseStrand = reverseStrand
        self.queryStartPosition = queryStartPosition
        self.queryEndPosition = queryEndPosition
        self.referenceStartPosition = referenceStartPosition
        self.referenceEndPosition = referenceEndPosition
        self.reverseStrand = reverseStrand
        self.confidence = confidence

    @property
    def positions(self):
        return [position for segment in self.segments for position in segment.positions]

    @property
    def alignedPairs(self) -> List[AlignedPair]:
        return [p for p in self.positions if isinstance(p, AlignedPair)]

    @property
    def notAlignedPositions(self) -> List[NotAlignedPosition]:
        return [p for p in self.positions if isinstance(p, NotAlignedPosition)]

    @property
    def cigarString(self):
        hitEnums = list(self.__getHitEnums())
        return "".join(self.__aggregateHitEnums(hitEnums))

    def __getHitEnums(self):
        pairs = list(self.__removeDuplicateQueryPositionsPreservingLastOne(self.alignedPairs))
        pairsIterator = iter(pairs)
        currentPair: AlignedPair = next(pairsIterator)
        previousQuery = currentPair.query.siteId
        for referenceIndex in range(pairs[0].reference.siteId,
                                    pairs[-1].reference.siteId + 1):
            queryIncrement = abs(currentPair.query.siteId - previousQuery)
            if queryIncrement > 1:
                for _ in range(1, queryIncrement):
                    yield HitEnum.INSERTION
                previousQuery = currentPair.query.siteId
            if currentPair.reference.siteId == referenceIndex:
                previousQuery = currentPair.query.siteId
                currentPair = next(pairsIterator, None)
                yield HitEnum.MATCH
            elif currentPair.reference.siteId > referenceIndex:
                yield HitEnum.DELETION

    @staticmethod
    def __removeDuplicateQueryPositionsPreservingLastOne(pairs: List[AlignedPair]):
        for _, ambiguousPairs in itertools.groupby(pairs, lambda pair: pair.query.siteId):
            *_, lastPair = ambiguousPairs
            yield lastPair

    @staticmethod
    def __aggregateHitEnums(hits: List[HitEnum]):
        hit = None
        count = 1
        previousHit: HitEnum = hits[0]
        for hit in hits[1:]:
            if hit == previousHit:
                count += 1
            else:
                yield AlignmentResultRow.__hitToString(count, previousHit)
                previousHit = hit
                count = 1
        if hit:
            yield AlignmentResultRow.__hitToString(count, hit)

    @staticmethod
    def __hitToString(count, hit):
        x = f"{count}{hit.value}"
        return x
