from __future__ import annotations


class Peak:
    null: Peak

    def __init__(self, position: int, height: float, leftBase: int = 0, rightBase: int = 0) -> None:
        self.position = position
        self.height = height
        self.leftProminenceBasePosition = leftBase
        self.rightProminenceBasePosition = rightBase


Peak.null = Peak(0, 0)
