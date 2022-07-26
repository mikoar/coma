from abc import ABC
from typing import NamedTuple, List


class XmapAlignmentPosition(NamedTuple):
    siteId: int


class XmapAlignedPair(NamedTuple):
    reference: XmapAlignmentPosition
    query: XmapAlignmentPosition


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
