from __future__ import annotations

import dataclasses
import itertools
from dataclasses import dataclass
from typing import List, Callable, Iterable

from src.alignment.aligned_pair import AlignedPair, HitEnum
from src.alignment.region_score_penalties import RegionScorePenalty
from src.alignment.region_scores import RegionScores


@dataclass
class AlignmentResults:
    referenceFilePath: str
    queryFilePath: str
    rows: List[AlignmentResultRow]


@dataclass
class AlignmentResultRow:
    alignedPairs: List[AlignedPair]
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
    def confidence(self):
        return 6.66

    @property
    def cigarString(self):
        hitEnums = list(self.__getHitEnums())
        return "".join(self.__aggregateHitEnums(hitEnums))

    def merge(self, other: AlignmentResultRow):
        pairs = self.alignedPairs + other.alignedPairs
        deduplicatedPairs = self.__deduplicatePairs(self.__deduplicatePairs(pairs, AlignedPair.querySelector),
                                                    AlignedPair.referenceSelector)
        return dataclasses.replace(self, alignedPairs=list(deduplicatedPairs))

    @staticmethod
    def __deduplicatePairs(pairs: Iterable[AlignedPair], key: Callable[[AlignedPair], int]):
        sortedPairs = sorted(pairs, key=key)
        for _, ambiguousPairs in itertools.groupby(sortedPairs, key):
            yield min(ambiguousPairs, key=AlignedPair.distanceSelector)

    def getRegionScores(self, penalties: List[RegionScorePenalty], perfectMatchScore: int = 10000):
        return RegionScores(list(self.__getRegionScoresGenerator(penalties, perfectMatchScore)))

    def __getHitEnums(self):
        pairs = list(self.__removeDuplicateQueryPositionsPreservingLastOne(self.alignedPairs))
        pairsIterator = iter(pairs)
        currentPair: AlignedPair = next(pairsIterator)
        previousQuery = currentPair.queryPositionIndex
        for referenceIndex in range(pairs[0].referencePositionIndex,
                                    pairs[-1].referencePositionIndex + 1):
            queryIncrement = abs(currentPair.queryPositionIndex - previousQuery)
            if queryIncrement > 1:
                for _ in range(1, queryIncrement):
                    yield HitEnum.INSERTION
                previousQuery = currentPair.queryPositionIndex
            if currentPair.referencePositionIndex == referenceIndex:
                previousQuery = currentPair.queryPositionIndex
                currentPair = next(pairsIterator, None)
                yield HitEnum.MATCH
            elif currentPair.referencePositionIndex > referenceIndex:
                yield HitEnum.DELETION

    @staticmethod
    def __removeDuplicateQueryPositionsPreservingLastOne(pairs: List[AlignedPair]):
        for _, ambiguousPairs in itertools.groupby(pairs, lambda pair: pair.queryPositionIndex):
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

    def __getRegionScoresGenerator(self, penalties: List[RegionScorePenalty], perfectMatchScore: int):
        previousPair: AlignedPair | None = None
        for pair in self.alignedPairs:
            yield perfectMatchScore - sum([p.getPenalty(previousPair, pair) for p in penalties])
            previousPair = pair
