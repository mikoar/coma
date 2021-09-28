from __future__ import annotations

from dataclasses import dataclass
from typing import List, Iterable

from src.alignment.aligned_pair import AlignedPair
from src.alignment.alignment_score import RegionScorePenalties


@dataclass
class AlignmentResult:
    referenceStartPosition: int
    referenceEndPosition: int
    alignedPairs: List[AlignedPair]

    def getRegionScores(self, penalties: RegionScorePenalties, perfectMatchScore: int = 10000) -> Iterable[int | float]:
        previousPair: AlignedPair | None = None
        for pair in self.alignedPairs:
            yield (perfectMatchScore
                   - penalties.getUnmatchedLabelPenalty(previousPair, pair)
                   - penalties.getDistancePenalty(previousPair, pair))
            previousPair = pair
