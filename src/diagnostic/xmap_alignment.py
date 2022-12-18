from __future__ import annotations

from abc import ABC
from dataclasses import dataclass, field
from typing import NamedTuple, List


class XmapAlignmentPosition(NamedTuple):
    siteId: int
    position: int


@dataclass(frozen=True)
class XmapAlignedPair:
    reference: XmapAlignmentPosition
    query: XmapAlignmentPosition

    @staticmethod
    def create(reference: str, query: str):
        return XmapAlignedPair(XmapAlignmentPosition(int(reference), 0), XmapAlignmentPosition(int(query), 0))

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId})"


@dataclass(frozen=True)
class XmapAlignedPairWithDistance(XmapAlignedPair):
    distance: int = field(compare=False)

    @staticmethod
    def calculateDistance(pair: XmapAlignedPair, firstPair: XmapAlignedPair | None, reverseStrand: bool):
        def queryDifference():
            return firstPair.query.position - pair.query.position if reverseStrand else pair.query.position - firstPair.query.position

        distance = queryDifference() - (pair.reference.position - firstPair.reference.position) if firstPair else 0
        return XmapAlignedPairWithDistance(pair.reference, pair.query, distance)

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {int(self.reference.position)}, {self.query.siteId}, " \
               f"{int(self.reference.position)}, {round(self.distance)})"


class XmapAlignment(ABC):
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
    alignedPairs: List[XmapAlignedPair]
    null: XmapAlignment

    @property
    def orientation(self):
        return "-" if self.reverseStrand else "+"


class _NullXmapAlignment(XmapAlignment):
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
    alignedPairs: List[XmapAlignedPair] = []

    @property
    def orientation(self):
        return "None"


XmapAlignment.null = _NullXmapAlignment()
