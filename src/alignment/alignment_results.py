from __future__ import annotations

from dataclasses import dataclass
from itertools import zip_longest
from typing import List

from src.alignment.aligned_pair import AlignedPair, HitEnum
from src.alignment.region_score_penalties import RegionScorePenalties
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
    def score(self):
        return 0.

    @property
    def cigarString(self):
        hits = [hit for pair, nextPair in zip_longest(self.alignedPairs, self.alignedPairs[1:]) for hit in
                pair.getHitEnums(nextPair)]
        return "".join(self.__hitEnumGenerator(hits))

    @staticmethod
    def __hitEnumGenerator(hits: List[HitEnum]):
        count = 1
        previousHit: HitEnum = hits[0]
        for hit in hits[1:]:
            if hit == previousHit:
                count += 1
            else:
                yield AlignmentResultRow.__hitToString(count, previousHit)
                previousHit = hit
                count = 1
        yield AlignmentResultRow.__hitToString(count, hit)

    @staticmethod
    def __hitToString(count, hit):
        x = f"{count}{hit.value}"
        return x

    def getRegionScores(self, penalties: RegionScorePenalties, perfectMatchScore: int = 10000):
        return RegionScores(list(self.__getRegionScoresGenerator(penalties, perfectMatchScore)))

    def __getRegionScoresGenerator(self, penalties: RegionScorePenalties, perfectMatchScore: int):
        previousPair: AlignedPair | None = None
        for pair in self.alignedPairs:
            yield (perfectMatchScore
                   - penalties.getUnmatchedLabelPenalty(previousPair, pair)
                   - penalties.getDistancePenalty(previousPair, pair))
            previousPair = pair
