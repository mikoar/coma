from cmap_reader import Alignment
from optical_map import CorrelationResult


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, result: CorrelationResult, reference: Alignment):
        return reference.refStartPosition <= result.peaks.max * self.resolution <= reference.refEndPosition
