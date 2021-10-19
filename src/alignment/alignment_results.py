from __future__ import annotations

from dataclasses import dataclass
from typing import List

from src.alignment.aligned_pair import AlignedPair
from src.alignment.region_score_penalties import RegionScorePenalties
from src.alignment.region_scores import RegionScores


@dataclass
class AlignmentResults:
    referenceFilePath: str
    queryFilePath: str
    rows: List[AlignmentResultRow]


@dataclass
class AlignmentResultRow:
    queryId: int
    referenceId: int
    referenceStartPosition: int
    referenceEndPosition: int
    reverseStrand: bool
    alignedPairs: List[AlignedPair]

    @property
    def score(self):
        return 0.

    def getRegionScores(self, penalties: RegionScorePenalties, perfectMatchScore: int = 10000):
        return RegionScores(list(self.__getRegionScoresGenerator(penalties, perfectMatchScore)))

    def __getRegionScoresGenerator(self, penalties: RegionScorePenalties, perfectMatchScore: int):
        previousPair: AlignedPair | None = None
        for pair in self.alignedPairs:
            yield (perfectMatchScore
                   - penalties.getUnmatchedLabelPenalty(previousPair, pair)
                   - penalties.getDistancePenalty(previousPair, pair))
            previousPair = pair
