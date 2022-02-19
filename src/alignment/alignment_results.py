from __future__ import annotations

import dataclasses
import itertools
from dataclasses import dataclass
from typing import List, Callable, Iterable

from src.alignment.alignment_position import AlignedPair, HitEnum, NotAlignedPosition, AlignmentPosition


@dataclass
class AlignmentResults:
    referenceFilePath: str
    queryFilePath: str
    rows: List[AlignmentResultRow]


@dataclass
class AlignmentResultRow:
    positions: List[AlignmentPosition]
    queryId: int = 1
    referenceId: int = 1
    queryStartPosition: int = 1
    queryEndPosition: int = 1
    referenceStartPosition: int = 1
    referenceEndPosition: int = 1
    queryLength: int = 1
    referenceLength: int = 1
    reverseStrand: bool = False

    @property
    def alignedPairs(self) -> List[AlignedPair]:
        return [p for p in self.positions if isinstance(p, AlignedPair)]

    @property
    def notAlignedPositions(self) -> List[NotAlignedPosition]:
        return [p for p in self.positions if isinstance(p, NotAlignedPosition)]

    @property
    def confidence(self):
        return 6.66

    @property
    def cigarString(self):
        hitEnums = list(self.__getHitEnums())
        return "".join(self.__aggregateHitEnums(hitEnums))

    def merge(self, other: AlignmentResultRow):
        pairs = self.alignedPairs + other.alignedPairs
        deduplicatedPairs = self.__deduplicatePairs(self.__deduplicatePairs(pairs, AlignedPair.querySiteIdSelector),
                                                    AlignedPair.referenceSiteIdSelector)
        return dataclasses.replace(self, positions=list(deduplicatedPairs))

    @staticmethod
    def __deduplicatePairs(pairs: Iterable[AlignedPair], key: Callable[[AlignedPair], int]):
        sortedPairs = sorted(pairs, key=key)
        for _, ambiguousPairs in itertools.groupby(sortedPairs, key):
            yield min(ambiguousPairs, key=AlignedPair.distanceSelector)

    def getScores(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int):
        return [p.getScore(perfectMatchScore, scoreMultiplier, unmatchedPenalty) for p in self.positions]

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
