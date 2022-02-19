from typing import List

from src.correlation.vectorise import vectorisePositions, blur


class SequenceGenerator:
    def __init__(self, resolution: int, blurRadius: int) -> None:
        self.resolution = resolution
        self.blurRadius = blurRadius

    def positionsToSequence(self, positions: List, start: int = 0, end: int = None):
        vector = list(vectorisePositions(positions, self.resolution, start, end))
        return blur(vector, self.blurRadius)
