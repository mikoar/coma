from __future__ import annotations

from dataclasses import dataclass
from typing import List

from src.alignment.aligned_pair import AlignedPair
from src.alignment.region_score_penalties import RegionScorePenalties
from src.alignment.region_scores import RegionScores


@dataclass
class AlignmentResult:
    referenceStartPosition: int
    referenceEndPosition: int
    alignedPairs: List[AlignedPair]

    def getRegionScores(self, penalties: RegionScorePenalties, perfectMatchScore: int = 10000):
        return RegionScores(list(self.__getRegionScoresGenerator(penalties, perfectMatchScore)))

    def __getRegionScoresGenerator(self, penalties, perfectMatchScore):
        previousPair: AlignedPair | None = None
        for pair in self.alignedPairs:
            yield (perfectMatchScore
                   - penalties.getUnmatchedLabelPenalty(previousPair, pair)
                   - penalties.getDistancePenalty(previousPair, pair))
            previousPair = pair
