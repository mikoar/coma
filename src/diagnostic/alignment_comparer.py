from __future__ import annotations

import statistics
from dataclasses import dataclass
from difflib import SequenceMatcher
from enum import Enum
from itertools import groupby
from typing import List, Dict, TextIO

import pandas as pd
from pandas import DataFrame

from src.diagnostic.benchmark_alignment import BenchmarkAlignment, BenchmarkAlignedPair


class AlignmentComparison:
    avgOverlappingAlignment1Coverage: float
    avgOverlappingAlignment2Coverage: float
    avgOverlappingIdentity: float
    overlapping: int
    nonOverlapping: int
    firstOnly: int
    secondOnly: int
    rows: List[AlignmentRowComparison]
    null: AlignmentComparison

    def __init__(self, avgAlignment1Coverage: float,
                 avgAlignment2Coverage: float,
                 avgIdentity: float,
                 overlapping: int,
                 nonOverlapping: int,
                 firstOnly: int,
                 secondOnly: int,
                 rows: List[AlignmentRowComparison]):
        self.avgOverlappingAlignment1Coverage = avgAlignment1Coverage
        self.avgOverlappingAlignment2Coverage = avgAlignment2Coverage
        self.avgOverlappingIdentity = avgIdentity
        self.overlapping = overlapping
        self.nonOverlapping = nonOverlapping
        self.firstOnly = firstOnly
        self.secondOnly = secondOnly
        self.rows = rows

    @staticmethod
    def create(rows: List[AlignmentRowComparison]):
        if not rows:
            return AlignmentComparison.null

        overlappingRows = [row for row in rows if row.overlapping]
        return AlignmentComparison(
            statistics.fmean(map(lambda row: row.alignment1Coverage, overlappingRows)) if overlappingRows else 0,
            statistics.fmean(map(lambda row: row.alignment2Coverage, overlappingRows)) if overlappingRows else 0,
            statistics.fmean(map(lambda row: row.identity, overlappingRows)) if overlappingRows else 0,
            len(overlappingRows),
            sum(1 for row in rows if row.type == AlignmentRowComparisonResultType.BOTH and not row.overlapping),
            sum(1 for row in rows if row.type == AlignmentRowComparisonResultType.FIRST_ONLY),
            sum(1 for row in rows if row.type == AlignmentRowComparisonResultType.SECOND_ONLY),
            rows)

    def write(self, file: TextIO, includePositions: bool):
        file.writelines([
            f"# AvgOverlappingAlignment1Coverage\t{self.avgOverlappingAlignment1Coverage}\n",
            f"# AvgOverlappingAlignment2Coverage\t{self.avgOverlappingAlignment2Coverage}\n",
            f"# AvgOverlappingIdentity\t{self.avgOverlappingIdentity}\n",
            f"# Overlapping\t{self.overlapping}\n",
            f"# NonOverlapping\t{self.nonOverlapping}\n",
            f"# FirstOnly\t{self.firstOnly}\n",
            f"# SecondOnly\t{self.secondOnly}\n"
        ])

        alignmentHeaderDescription = "(referenceID, referencePosition, queryID, queryPosition, distance)" \
            if includePositions else "(referenceID, queryID)"

        headers = [
            "QryContigID",
            "RefContigID",
            "Type",
            "Identity",
            "Alignment1Coverage",
            "Alignment2Coverage",
            "Orientation",
            f"Alignment1Diff {alignmentHeaderDescription}",
            f"Alignment2Diff {alignmentHeaderDescription}",
            f"Alignment1 {alignmentHeaderDescription}",
            f"Alignment2 {alignmentHeaderDescription}"
        ]
        file.write("\t".join([header for header in ["#"] + headers]) + "\n")

        data = [[
            row.queryId,
            row.referenceId,
            row.type.name,
            "{:.3f}".format(row.identity),
            "{:.3f}".format(row.alignment1Coverage),
            "{:.3f}".format(row.alignment2Coverage),
            row.orientation,
            "".join([r.toString(includePositions) for r in row.alignment1ExclusivePairs]) if row.overlapping else "",
            "".join([r.toString(includePositions) for r in row.alignment2ExclusivePairs]) if row.overlapping else "",
            "".join([r.toString(includePositions) for r in row.alignment1.alignedPairs]),
            "".join([r.toString(includePositions) for r in row.alignment2.alignedPairs])
        ] for row in self.rows]

        dataFrame = DataFrame(data, columns=headers, index=pd.RangeIndex(start=1, stop=len(self.rows) + 1))
        dataFrame.to_csv(file, sep='\t', header=False, mode="a", lineterminator="\n")


class _NullAlignmentComparison(AlignmentComparison):
    def __init__(self):
        super().__init__(0., 0., 0., 0, 0, 0, 0, [])

    def write(self, file: TextIO, includePositions: bool):
        return


AlignmentComparison.null = _NullAlignmentComparison()


class AlignmentRowComparisonResultType(Enum):
    BOTH = 1,
    FIRST_ONLY = 2,
    SECOND_ONLY = 3


@dataclass
class AlignmentRowComparison:
    type: AlignmentRowComparisonResultType
    alignment1: BenchmarkAlignment
    alignment2: BenchmarkAlignment
    alignment1ExclusivePairs: List[BenchmarkAlignedPair]
    alignment2ExclusivePairs: List[BenchmarkAlignedPair]
    alignment1Coverage: float
    alignment2Coverage: float
    identity: float

    @property
    def queryId(self):
        return self.alignment1.queryId or self.alignment2.queryId

    @property
    def referenceId(self):
        return self.alignment1.referenceId or self.alignment2.referenceId

    @property
    def orientation(self):
        return self.alignment1.orientation \
            if self.alignment1.orientation == self.alignment2.orientation \
            else f"{self.alignment1.orientation}/{self.alignment2.orientation}"

    @property
    def overlapping(self):
        return self.identity > 0.

    @staticmethod
    def alignment1Only(alignment1: BenchmarkAlignment):
        return AlignmentRowComparison(
            AlignmentRowComparisonResultType.FIRST_ONLY, alignment1, BenchmarkAlignment.null, [], [], 0., 0., 0.)

    @staticmethod
    def alignment2Only(alignment2: BenchmarkAlignment):
        return AlignmentRowComparison(
            AlignmentRowComparisonResultType.SECOND_ONLY, BenchmarkAlignment.null, alignment2, [], [], 0., 0., 0.)


class AlignmentComparer:
    def __init__(self, rowComparer: AlignmentRowComparer):
        self.__rowComparer = rowComparer

    def compare(self, alignments1: List[BenchmarkAlignment], alignments2: List[BenchmarkAlignment]):
        alignments1Dict = self.__toDict(alignments1)
        alignments2Dict = self.__toDict(alignments2)
        comparedRows = [self.__rowComparer.compare(a1, alignments2Dict[key]) for key, a1 in alignments1Dict.items() if
                        key in alignments2Dict]
        alignments1OnlyRows = [AlignmentRowComparison.alignment1Only(a1) for a1
                               in self.__getNotMatchingAlignments(alignments1Dict, alignments2Dict)]
        alignments2OnlyRows = [AlignmentRowComparison.alignment2Only(a2) for a2
                               in self.__getNotMatchingAlignments(alignments2Dict, alignments1Dict)]
        return AlignmentComparison.create(comparedRows + alignments1OnlyRows + alignments2OnlyRows)

    @staticmethod
    def __toDict(alignments: List[BenchmarkAlignment]) -> Dict[(int, int), BenchmarkAlignment]:
        return {(a.queryId, a.referenceId): a for a in sorted(alignments, key=lambda a: (a.referenceId, a.queryId))}

    @staticmethod
    def __getNotMatchingAlignments(source: Dict[(int, int), BenchmarkAlignment], target: Dict[(int, int), BenchmarkAlignment]):
        return [a1 for key, a1 in source.items() if key not in target]


class AlignmentRowComparer:
    def __init__(self, combineMultipleQuerySources: bool):
        self.combineMultipleQuerySources = combineMultipleQuerySources

    def compare(self, alignment1: BenchmarkAlignment, alignment2: BenchmarkAlignment):
        alignment1Pairs = self.__combineMultipleQuerySources(alignment1.alignedPairs, alignment2.alignedPairs)
        alignment2Pairs = self.__combineMultipleQuerySources(alignment2.alignedPairs, alignment1.alignedPairs)
        difference1 = self.__getDifference(alignment1Pairs, alignment2Pairs)
        coverage1 = self.__getCoverage(alignment1Pairs, difference1)
        difference2 = self.__getDifference(alignment2Pairs, alignment1Pairs)
        coverage2 = self.__getCoverage(alignment2Pairs, difference2)

        ratio = self.__getIdentityRatio(alignment1Pairs, alignment2Pairs)
        return AlignmentRowComparison(
            AlignmentRowComparisonResultType.BOTH, alignment1, alignment2, difference1, difference2, coverage1, coverage2, ratio)

    def __combineMultipleQuerySources(self, pairs: List[BenchmarkAlignedPair], otherPairs: List[BenchmarkAlignedPair]):
        if not self.combineMultipleQuerySources:
            return pairs

        alignmentsPerQuery = [self.__removeExtraAlignmentsWithSameQueryIfOneOfThemIsInOtherList(list(alignmentsWithSameQuery), otherPairs)
                              for _, alignmentsWithSameQuery in groupby(pairs, BenchmarkAlignedPair.querySiteIdSelector)]
        return [a for alignments in alignmentsPerQuery for a in alignments]

    @staticmethod
    def __getIdentityRatio(alignment1Pairs: List[BenchmarkAlignedPair], alignment2Pairs: List[BenchmarkAlignedPair]):
        matcher = SequenceMatcher(None, alignment1Pairs, alignment2Pairs)
        ratio = matcher.ratio()
        return ratio

    @staticmethod
    def __getDifference(pairs: List[BenchmarkAlignedPair], otherPairs: List[BenchmarkAlignedPair]):
        return sorted(set(pairs).difference(set(otherPairs)), key=BenchmarkAlignedPair.referenceSiteIdSelector)

    @staticmethod
    def __removeExtraAlignmentsWithSameQueryIfOneOfThemIsInOtherList(
            alignmentsWithSameQuery: List[BenchmarkAlignedPair], otherPairs: List[BenchmarkAlignedPair]):
        return [a for a in alignmentsWithSameQuery if a in otherPairs] or alignmentsWithSameQuery

    @staticmethod
    def __getCoverage(pairs: List[BenchmarkAlignedPair], difference: List[BenchmarkAlignedPair]):
        pairsLength = len(pairs)
        return (pairsLength - len(difference)) / pairsLength if pairsLength > 0 else 1.
