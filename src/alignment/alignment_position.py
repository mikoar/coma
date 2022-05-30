from __future__ import annotations

import itertools
from abc import ABC, abstractmethod
from typing import Tuple, Iterable, Callable

from src.correlation.optical_map import PositionWithSiteId


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
        if unmatchedPenalty > 0:
            raise ValueError("penalty should be negative")
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
    null: AlignedPair

    def __init__(self, reference: PositionWithSiteId, query: PositionWithSiteId, queryShift: int = 0, source: int = 0):
        self.reference = reference
        self.query = query
        self.queryShift = queryShift
        self.source = source

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

    @staticmethod
    def deduplicate(pairs: Iterable[AlignedPair]):
        return AlignedPair.__deduplicateByKey(
            AlignedPair.__deduplicateByKey(pairs, AlignedPair.querySiteIdSelector),
            AlignedPair.referenceSiteIdSelector)

    @staticmethod
    def __deduplicateByKey(pairs: Iterable[AlignedPair], key: Callable[[AlignedPair], int]):
        sortedPairs = sorted(pairs, key=key)
        for _, ambiguousPairs in itertools.groupby(sortedPairs, key):
            yield min(ambiguousPairs, key=AlignedPair.distanceSelector)

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

    def lessOnBothSequences(self, other: AlignedPair):
        return self.query < other.query and self.reference < other.reference

    def lessOrEqualOnAnySequence(self, other: AlignedPair):
        return self.query < other.query or self.reference < other.reference \
               or self.query == other.query or self.reference == other.reference

    def __repr__(self) -> str:
        return f"({self.reference.siteId}, {self.query.siteId}), distance:{self.queryShift:.2f}, source:{self.source}"

    def __eq__(self, other: AlignedPair | Tuple[int, int] | Tuple[int, int, int]) -> bool:
        if isinstance(other, AlignedPair):
            return self.query == other.query and self.reference == other.reference
        return self.reference.siteId == other[0] and self.query.siteId == other[1] and (
                len(other) == 2 or self.queryShift == other[2])


AlignedPair.null = AlignedPair(PositionWithSiteId(0, 0), PositionWithSiteId(0, 0))


class ScoredAlignmentPosition(AlignmentPosition, ABC):
    score: float


class ScoredAlignedPair(AlignedPair, ScoredAlignmentPosition):
    def __init__(self, pair: AlignedPair, score: float):
        super().__init__(pair.reference, pair.query, pair.queryShift, pair.source)
        self.score = score

    def __repr__(self) -> str:
        return f"{AlignedPair.__repr__(self)}, score:{self.score:.2f}"


class ScoredNotAlignedPosition(NotAlignedPosition, ScoredAlignmentPosition):
    @property
    def absolutePosition(self) -> int:
        return self.__position.absolutePosition

    def __init__(self, position: NotAlignedPosition, score: float):
        self.__position = position
        self.score = score

    def __repr__(self) -> str:
        return f"{self.__position}, score:{self.score:.2f}"

    def __eq__(self, other: NotAlignedPosition):
        return self.__position == other
