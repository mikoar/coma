from __future__ import annotations

from dataclasses import dataclass
from difflib import SequenceMatcher
from enum import Enum
from typing import Tuple, List, Dict, TextIO

import pandas as pd
from pandas import DataFrame

from src.diagnostic.xmap_alignment import XmapAlignment


@dataclass
class AlignmentComparison:
    rows: List[AlignmentRowComparison]

    def write(self, file: TextIO):
        if not self.rows:
            return

        headers = [
            "QryContigID",
            "RefContigID",
            "Type",
            "Identity",
            "Alignment1Coverage",
            "Alignment2Coverage",
            "Alignment1",
            "Alignment2"
        ]
        file.write("\t".join([header for header in headers]) + "\n")

        data = [[
            row.queryId,
            row.referenceId,
            row.type.name,
            "{:.3f}".format(row.identity),
            "{:.3f}".format(row.alignment1Coverage),
            "{:.3f}".format(row.alignment2Coverage),
            "".join(str(row.pairs1)),
            "".join(str(row.pairs2))
        ] for row in self.rows]

        dataFrame = DataFrame(data, columns=headers, index=pd.RangeIndex(start=1, stop=len(self.rows) + 1))
        dataFrame.to_csv(file, sep='\t', header=False, mode="a", line_terminator="\n")


class AlignmentRowComparisonResultType(Enum):
    BOTH = 1,
    FIRST_ONLY = 2,
    SECOND_ONLY = 3


@dataclass
class AlignmentRowComparison:
    queryId: int
    referenceId: int
    type: AlignmentRowComparisonResultType
    pairs1: List[Tuple[int, int]]
    alignment1Coverage: float
    pairs2: List[Tuple[int, int]]
    alignment2Coverage: float
    identity: float

    @staticmethod
    def alignment1Only(queryId: int, referenceId: int, pairs1: List[Tuple[int, int]]):
        return AlignmentRowComparison(queryId, referenceId, AlignmentRowComparisonResultType.FIRST_ONLY, pairs1, 0., [],
                                      0., 0.)

    @staticmethod
    def alignment2Only(queryId: int, referenceId: int, pairs2: List[Tuple[int, int]]):
        return AlignmentRowComparison(queryId, referenceId, AlignmentRowComparisonResultType.SECOND_ONLY, [], 0.,
                                      pairs2, 0., 0.)


class AlignmentComparer:
    def __init__(self, rowComparer: AlignmentRowComparer):
        self.__rowComparer = rowComparer

    def compare(self, alignments1: List[XmapAlignment], alignments2: List[XmapAlignment]):
        alignments1Dict = self.__toDict(alignments1)
        alignments2Dict = self.__toDict(alignments2)
        comparedRows = [self.__rowComparer.compare(a1, alignments2Dict[key]) for key, a1 in alignments1Dict.items() if
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
        alignment1Coverage = self.__getCoverage(pairs1, pairs2)
        alignment2Coverage = self.__getCoverage(pairs2, pairs1)
        matcher = SequenceMatcher(None, pairs1, pairs2)
        ratio = matcher.ratio()
        return AlignmentRowComparison(referenceAlignmentRow.queryId, referenceAlignmentRow.referenceId,
                                      AlignmentRowComparisonResultType.BOTH, pairs1, alignment1Coverage, pairs2,
                                      alignment2Coverage, ratio)

    @staticmethod
    def __getCoverage(pairs: List[Tuple[int, int]], otherPairs: List[Tuple[int, int]]):
        pairsLength = len(pairs)
        return (pairsLength - len(set(pairs).difference(set(otherPairs)))) / pairsLength if pairsLength > 0 else 1.
