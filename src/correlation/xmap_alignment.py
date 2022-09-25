from abc import ABC
from typing import NamedTuple, List


class XmapAlignmentPosition(NamedTuple):
    siteId: int
    position: int


class XmapAlignedPair(NamedTuple):
    reference: XmapAlignmentPosition
    query: XmapAlignmentPosition

    @staticmethod
    def create(reference: str, query: str):
        return XmapAlignedPair(XmapAlignmentPosition(int(reference), 0), XmapAlignmentPosition(int(query), 0))

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId})"


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

    @property
    def orientation(self):
        return "-" if self.reverseStrand else "+"
