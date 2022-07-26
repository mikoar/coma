from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple, List

from src.correlation.xmap_alignment import XmapAlignment


@dataclass
class AlignmentComparisonResult:
    pairs1: List[Tuple[int, int]]
    query1Coverage: float
    pairs2: List[Tuple[int, int]]
    query2Coverage: float
    queryId: int
    referenceId: int


class AlignmentComparer:
    def compare(self, referenceAlignment: XmapAlignment, actualAlignment: XmapAlignment):
        pairs1 = list(
            map(lambda pair: (int(pair.reference.siteId), int(pair.query.siteId)), referenceAlignment.alignedPairs))
        pairs2 = list(
            map(lambda pair: (int(pair.reference.siteId), int(pair.query.siteId)), actualAlignment.alignedPairs))
        query1Coverage = self.__getCoverage(pairs1, pairs2)
        query2Coverage = self.__getCoverage(pairs2, pairs1)
        return AlignmentComparisonResult(pairs1, query1Coverage, pairs2, query2Coverage, referenceAlignment.queryId,
                                         referenceAlignment.referenceId)

    @staticmethod
    def __getCoverage(pairs: List[Tuple[int, int]], otherPairs: List[Tuple[int, int]]):
        pairsLength = len(pairs)
        return (pairsLength - len(set(pairs).difference(set(otherPairs)))) / pairsLength if pairsLength > 0 else 1.
