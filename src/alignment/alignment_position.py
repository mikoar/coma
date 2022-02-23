from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import Tuple, Final

from src.correlation.optical_map import PositionWithSiteId


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"


class AlignmentPosition(ABC):
    @abstractmethod
    def getScoredPosition(self, perfectMatchScore: int, scoreMultiplier: int,
                          unmatchedPenalty: int) -> ScoredAlignmentPosition:
        pass

    @property
    @abstractmethod
    def absolutePosition(self) -> int:
        pass

    def __lt__(self, other: AlignmentPosition):
        return self.absolutePosition < other.absolutePosition


class NotAlignedPosition(AlignmentPosition, ABC):
    def getScoredPosition(self, perfectMatchScore: int, scoreMultiplier: int,
                          unmatchedPenalty: int) -> ScoredAlignmentPosition:
        return ScoredNotAlignedPosition(self, unmatchedPenalty)


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

    def getScoredPosition(self, perfectMatchScore: int, scoreMultiplier: int,
                          unmatchedPenalty: int) -> ScoredAlignmentPosition:
        score = scoreMultiplier * (perfectMatchScore - self.distance)
        return ScoredAlignedPair(self, score)

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId}), d:{self.queryShift:.2f}"

    def __eq__(self, other: Tuple[int, int] | Tuple[int, int, int]) -> bool:
        return self.reference.siteId == other[0] and self.query.siteId == other[1] and (
                len(other) == 2 or self.queryShift == other[2])


nullAlignedPair: Final[AlignedPair] = AlignedPair(PositionWithSiteId(0, 0), PositionWithSiteId(0, 0))


class ScoredAlignmentPosition(AlignmentPosition, ABC):
    score: float


class ScoredAlignedPair(AlignedPair, ScoredAlignmentPosition):
    def __init__(self, pair: AlignedPair, score: float):
        super().__init__(pair.reference, pair.query, pair.queryShift)
        self.score = score

    def __repr__(self) -> str:
        return f"{AlignedPair.__repr__(self)}, s:{self.score:.2f}"


class ScoredNotAlignedPosition(NotAlignedPosition, ScoredAlignmentPosition):
    @property
    def absolutePosition(self) -> int:
        return self.__position.absolutePosition

    def __init__(self, position: NotAlignedPosition, score: float):
        self.__position = position
        self.score = score

    def __repr__(self) -> str:
        return f"{self.__position}, s:{self.score:.2f}"
