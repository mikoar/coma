from __future__ import annotations

from enum import Enum
from typing import NamedTuple, Tuple


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"


class NotAlignedPosition:
    def __init__(self, queryPositionIndex: int):
        self.queryPositionIndex = queryPositionIndex

    def __repr__(self) -> str:
        return str(self.queryPositionIndex)

    def __eq__(self, other: NotAlignedPosition | int) -> bool:
        return self.queryPositionIndex == other if isinstance(other, int) \
            else self.queryPositionIndex == other.queryPositionIndex


class AlignedPair(NamedTuple):
    referencePositionIndex: int
    queryPositionIndex: int
    queryShift: int = 0

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

    def getFalseNegativesCount(self, previousPair: AlignedPair | None):
        return 0 if not previousPair else self.referencePositionIndex - previousPair.referencePositionIndex - 1

    def getFalsePositivesCount(self, previousPair: AlignedPair | None):
        return 0 if not previousPair else abs(self.queryPositionIndex - previousPair.queryPositionIndex) - 1

    def __str__(self) -> str:
        return f"({self.referencePositionIndex}, {self.queryPositionIndex})"

    def __repr__(self) -> str:
        return f"({self.referencePositionIndex}, {self.queryPositionIndex}, {self.queryShift:.2f})"

    def __eq__(self, other: Tuple[int, int] | Tuple[int, int, int]) -> bool:
        return self.referencePositionIndex == other[0] and self.queryPositionIndex == other[1] and (
                len(other) == 2 or self.queryShift == other[2])


nullAlignedPair = AlignedPair(0, 0, 0)
