from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import Tuple

from src.correlation.optical_map import PositionWithSiteId


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"


class AlignmentPosition(ABC):
    @abstractmethod
    def getScore(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int) -> float:
        pass

    @property
    @abstractmethod
    def absolutePosition(self) -> int:
        pass

    def __lt__(self, other: AlignmentPosition):
        return self.absolutePosition < other.absolutePosition


class NotAlignedPosition(AlignmentPosition, ABC):
    def getScore(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int) -> float:
        return unmatchedPenalty


class NotAlignedQueryPosition(NotAlignedPosition):
    def __init__(self, query: PositionWithSiteId, referenceStart: int):
        self.query = query
        self.referenceStart = referenceStart

    @property
    def absolutePosition(self) -> int:
        return self.query.position + self.referenceStart

    def __repr__(self) -> str:
        return f"(-, {self.query.siteId})"

    def __eq__(self, other: NotAlignedQueryPosition | Tuple[None, int]) -> bool:
        return other[0] is None and self.query.siteId == other[1] if len(other) == 2 \
            else isinstance(other, NotAlignedQueryPosition) \
                 and self.query.siteId == other.query.siteId


class NotAlignedReferencePosition(NotAlignedPosition):
    def __init__(self, reference: PositionWithSiteId):
        self.reference = reference

    @property
    def absolutePosition(self) -> int:
        return self.reference.position

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, -)"

    def __eq__(self, other: NotAlignedReferencePosition | Tuple[int, None]) -> bool:
        return self.reference.siteId == other[0] and other[1] is None if len(other) == 2 \
            else isinstance(other, NotAlignedReferencePosition) \
                 and self.reference.siteId == other.reference.siteId


class AlignedPair(AlignmentPosition):
    def __init__(self, reference: PositionWithSiteId, query: PositionWithSiteId, queryShift: int = 0):
        self.reference = reference
        self.query = query
        self.queryShift = queryShift

    @staticmethod
    def distanceSelector(pair: AlignedPair):
        return pair.distance

    @staticmethod
    def queryShiftSelector(pair: AlignedPair):
        return pair.queryShift

    @staticmethod
    def referenceSiteIdSelector(pair: AlignedPair):
        return pair.reference.siteId

    @staticmethod
    def querySiteIdSelector(pair: AlignedPair):
        return pair.query.siteId

    @property
    def distance(self):
        return abs(self.queryShift)

    @property
    def absolutePosition(self) -> int:
        return self.reference.position

    def getScore(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int) -> float:
        return scoreMultiplier * (perfectMatchScore - self.distance)

    def __str__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId})"

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId}, {self.queryShift:.2f})"

    def __eq__(self, other: Tuple[int, int] | Tuple[int, int, int]) -> bool:
        return self.reference.siteId == other[0] and self.query.siteId == other[1] and (
                len(other) == 2 or self.queryShift == other[2])


nullAlignedPair = AlignedPair(PositionWithSiteId(0, 0), PositionWithSiteId(0, 0))
