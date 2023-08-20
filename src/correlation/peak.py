from __future__ import annotations


class Peak:
    null: Peak

    def __init__(self, position: int, height: float, leftBase: int = 0, rightBase: int = 0, score: float = 0.) -> None:
        self.position = position
        self.height = height
        self.leftProminenceBasePosition = leftBase
        self.rightProminenceBasePosition = rightBase
        self.score = score

    @property
    def width(self):
        return self.rightProminenceBasePosition - self.leftProminenceBasePosition

    def __eq__(self, other):
        return isinstance(other, Peak) and self.position == other.position \
               and self.height == other.height \
               and self.leftProminenceBasePosition == other.leftProminenceBasePosition \
               and self.rightProminenceBasePosition == other.rightProminenceBasePosition


Peak.null = Peak(0, 0)
