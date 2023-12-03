from __future__ import annotations

from abc import ABC
from dataclasses import dataclass, field
from typing import NamedTuple, List


class BenchmarkAlignmentPosition(NamedTuple):
    siteId: int
    position: int


@dataclass(frozen=True)
class BenchmarkAlignedPair:
    reference: BenchmarkAlignmentPosition
    query: BenchmarkAlignmentPosition

    @staticmethod
    def create(reference: str, query: str):
        return BenchmarkAlignedPair(BenchmarkAlignmentPosition(int(reference), 0), BenchmarkAlignmentPosition(int(query), 0))

    @staticmethod
    def referenceSiteIdSelector(pair: BenchmarkAlignedPair):
        return pair.reference.siteId

    @staticmethod
    def querySiteIdSelector(pair: BenchmarkAlignedPair):
        return pair.query.siteId

    def toString(self, includePositions: bool):
        return self.__repr__()

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId})"


@dataclass(frozen=True)
class BenchmarkAlignedPairWithDistance(BenchmarkAlignedPair):
    distance: int = field(compare=False)

    @staticmethod
    def calculateDistance(pair: BenchmarkAlignedPair, firstPair: BenchmarkAlignedPair | None, reverseStrand: bool):
        def queryDifference():
            return firstPair.query.position - pair.query.position if reverseStrand else pair.query.position - firstPair.query.position

        distance = queryDifference() - (pair.reference.position - firstPair.reference.position) if firstPair else 0
        return BenchmarkAlignedPairWithDistance(pair.reference, pair.query, distance)

    def toString(self, includePositions: bool):
        return self.__repr__() if includePositions else super().__repr__()

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {int(self.reference.position)}, {self.query.siteId}, " \
               f"{int(self.reference.position)}, {round(self.distance)})"


class BenchmarkAlignment(ABC):
    queryId: int
    referenceId: int
    queryStartPosition: int
    queryEndPosition: int
    referenceStartPosition: int
    referenceEndPosition: int
    reverseStrand: bool
    confidence: float
    cigarString: str
    queryLength: int
    referenceLength: int
    alignedPairs: List[BenchmarkAlignedPair]
    null: BenchmarkAlignment

    @property
    def orientation(self):
        return "-" if self.reverseStrand else "+"


class _NullBenchmarkAlignment(BenchmarkAlignment):
    queryId: int = 0
    referenceId: int = 0
    queryStartPosition: int = 0
    queryEndPosition: int = 0
    referenceStartPosition: int = 0
    referenceEndPosition: int = 0
    reverseStrand: bool = False
    confidence: float = 0.
    cigarString: str = ""
    queryLength: int = 0
    referenceLength: int = 0
    alignedPairs: List[BenchmarkAlignedPair] = []

    @property
    def orientation(self):
        return "None"


BenchmarkAlignment.null = _NullBenchmarkAlignment()
