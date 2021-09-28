from __future__ import annotations

from typing import NamedTuple


class AlignedPair(NamedTuple):
    referencePositionIndex: int
    queryPositionIndex: int
    queryShift: int = 0

    @property
    def distance(self):
        return abs(self.queryShift)

    def getFalseNegativesCount(self, previousPair: AlignedPair | None):
        return 0 if not previousPair else self.referencePositionIndex - previousPair.referencePositionIndex - 1

    def getFalsePositivesCount(self, previousPair: AlignedPair | None):
        return 0 if not previousPair else abs(self.queryPositionIndex - previousPair.queryPositionIndex) - 1

    def __repr__(self) -> str:
        return f"({self.referencePositionIndex}, {self.queryPositionIndex})"

    def __eq__(self, other) -> bool:
        return self.referencePositionIndex == other[0] and self.queryPositionIndex == other[1]
