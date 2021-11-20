from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Tuple, List

from src.alignment.aligned_pair import AlignedPair
from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.bionano_alignment import BionanoAlignment


@dataclass
class AlignmentComparisonResult:
    pairs1: List[Tuple[int, int]]
    query1Coverage: float
    pairs2: List[Tuple[int, int]]
    query2Coverage: float
    queryId: int
    referenceId: int


class AlignmentComparer:
    def compare(self, referenceAlignment: BionanoAlignment, actualAlignment: AlignmentResultRow):
        pairs1 = self.__getPairTuples(referenceAlignment)
        pairs2 = self.__getPairTuples(actualAlignment)
        query1Coverage = self.__getCoverage(pairs1, pairs2)
        query2Coverage = self.__getCoverage(pairs2, pairs1)
        return AlignmentComparisonResult(pairs1, query1Coverage, pairs2, query2Coverage, referenceAlignment.queryId,
                                         referenceAlignment.referenceId)

    @staticmethod
    def __getPairTuples(alignment: BionanoAlignment | AlignmentResultRow):
        queryPositionSelector: Callable[[AlignedPair], Tuple[int, int]] = lambda pair: (
            int(pair.referencePositionIndex), int(pair.queryPositionIndex))
        return list(map(queryPositionSelector, alignment.alignedPairs))

    @staticmethod
    def __getCoverage(pairs: List[Tuple[int, int]], otherPairs: List[Tuple[int, int]]):
        pairsLength = len(pairs)
        return (pairsLength - len(set(pairs).difference(set(otherPairs)))) / pairsLength if pairsLength > 0 else 1.
