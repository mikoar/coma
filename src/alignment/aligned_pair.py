from __future__ import annotations

from abc import ABC, abstractmethod
from enum import Enum
from typing import Tuple


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"


class AlignmentPosition(ABC):
    @abstractmethod
    def getScore(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int) -> float:
        pass


class NotAlignedPosition(AlignmentPosition, ABC):
    def getScore(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int) -> float:
        return unmatchedPenalty


class NotAlignedQueryPosition(NotAlignedPosition):
    def __init__(self, queryPositionIndex: int):
        self.queryPositionIndex = queryPositionIndex

    def __repr__(self) -> str:
        return f"(-, {self.queryPositionIndex})"

    def __eq__(self, other: NotAlignedQueryPosition | Tuple[None, int]) -> bool:
        return other[0] is None and self.queryPositionIndex == other[1] if len(other) == 2 \
            else isinstance(other, NotAlignedQueryPosition) \
                 and self.queryPositionIndex == other.queryPositionIndex


class NotAlignedReferencePosition(NotAlignedPosition):
    def __init__(self, referencePositionIndex: int):
        self.referencePositionIndex = referencePositionIndex

    def __repr__(self) -> str:
        return f"({self.referencePositionIndex}, -)"

    def __eq__(self, other: NotAlignedReferencePosition | Tuple[int, None]) -> bool:
        return self.referencePositionIndex == other[0] and other[1] is None if len(other) == 2 \
            else isinstance(other, NotAlignedReferencePosition) \
                 and self.referencePositionIndex == other.referencePositionIndex


class AlignedPair(AlignmentPosition):
    def __init__(self, referencePositionIndex: int, queryPositionIndex: int, queryShift: int = 0):
        self.referencePositionIndex = referencePositionIndex
        self.queryPositionIndex = queryPositionIndex
        self.queryShift = queryShift

    @staticmethod
    def distanceSelector(pair: AlignedPair):
        return pair.distance

    @staticmethod
    def queryShiftSelector(pair: AlignedPair):
        return pair.queryShift

    @staticmethod
    def referenceSelector(pair: AlignedPair):
        return pair.referencePositionIndex

    @staticmethod
    def querySelector(pair: AlignedPair):
        return pair.queryPositionIndex

    @property
    def distance(self):
        return abs(self.queryShift)

    def getScore(self, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int) -> float:
        return scoreMultiplier * (perfectMatchScore - self.distance)

    def __str__(self) -> str:
        return f"({self.referencePositionIndex}, {self.queryPositionIndex})"

    def __repr__(self) -> str:
        return f"({self.referencePositionIndex}, {self.queryPositionIndex}, {self.queryShift:.2f})"

    def __eq__(self, other: Tuple[int, int] | Tuple[int, int, int]) -> bool:
        return self.referencePositionIndex == other[0] and self.queryPositionIndex == other[1] and (
                len(other) == 2 or self.queryShift == other[2])


nullAlignedPair = AlignedPair(0, 0)
