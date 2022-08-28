from __future__ import annotations

from dataclasses import dataclass
from difflib import SequenceMatcher
from typing import Tuple, List, Dict

from src.correlation.xmap_alignment import XmapAlignment


@dataclass
class AlignmentComparison:
    rows: List[AlignmentRowComparison]


@dataclass
class AlignmentRowComparison:
    queryId: int
    referenceId: int
    pairs1: List[Tuple[int, int]]
    query1Coverage: float
    pairs2: List[Tuple[int, int]]
    query2Coverage: float
    identity: float

    @staticmethod
    def alignment1Only(queryId: int, referenceId: int, pairs1: List[Tuple[int, int]]):
        return AlignmentRowComparison(queryId, referenceId, pairs1, 0., [], 0., 0.)

    @staticmethod
    def alignment2Only(queryId: int, referenceId: int, pairs2: List[Tuple[int, int]]):
        return AlignmentRowComparison(queryId, referenceId, [], 0., pairs2, 0., 0.)


class AlignmentComparer:
    def __init__(self, rowComparer: AlignmentRowComparer):
        self.rowComparer = rowComparer

    def compare(self, alignments1: List[XmapAlignment], alignments2: List[XmapAlignment]):
        alignments1Dict = self.__toDict(alignments1)
        alignments2Dict = self.__toDict(alignments2)
        comparedRows = [self.rowComparer.compare(a1, alignments2Dict[key]) for key, a1 in alignments1Dict.items() if
                        key in alignments2Dict]
        alignments1OnlyRows = [AlignmentRowComparison.alignment1Only(a1.queryId, a1.referenceId, a1.alignedPairs) for a1
                               in self.__getNotMatchingAlignments(alignments1Dict, alignments2Dict)]
        alignments2OnlyRows = [AlignmentRowComparison.alignment2Only(a2.queryId, a2.referenceId, a2.alignedPairs) for a2
                               in self.__getNotMatchingAlignments(alignments2Dict, alignments1Dict)]
        return AlignmentComparison(comparedRows + alignments1OnlyRows + alignments2OnlyRows)

    @staticmethod
    def __toDict(alignments: List[XmapAlignment]) -> Dict[(int, int), XmapAlignment]:
        return {(a.queryId, a.referenceId): a for a in sorted(alignments, key=lambda a: (a.referenceId, a.queryId))}

    @staticmethod
    def __getNotMatchingAlignments(source: Dict[(int, int), XmapAlignment], target: Dict[(int, int), XmapAlignment]):
        return [a1 for key, a1 in source.items() if key not in target]


class AlignmentRowComparer:
    def compare(self, referenceAlignmentRow: XmapAlignment, actualAlignmentRow: XmapAlignment):
        pairs1 = list(
            map(lambda pair: (int(pair.reference.siteId), int(pair.query.siteId)), referenceAlignmentRow.alignedPairs))
        pairs2 = list(
            map(lambda pair: (int(pair.reference.siteId), int(pair.query.siteId)), actualAlignmentRow.alignedPairs))
        query1Coverage = self.__getCoverage(pairs1, pairs2)
        query2Coverage = self.__getCoverage(pairs2, pairs1)
        matcher = SequenceMatcher(None, pairs1, pairs2)
        ratio = matcher.ratio()
        return AlignmentRowComparison(referenceAlignmentRow.queryId, referenceAlignmentRow.referenceId, pairs1,
                                      query1Coverage, pairs2, query2Coverage, ratio)

    @staticmethod
    def __getCoverage(pairs: List[Tuple[int, int]], otherPairs: List[Tuple[int, int]]):
        pairsLength = len(pairs)
        return (pairsLength - len(set(pairs).difference(set(otherPairs)))) / pairsLength if pairsLength > 0 else 1.
