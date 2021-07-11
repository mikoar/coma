from cmap_reader import Alignment
from optical_map import CorrelationResult, Peaks


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, result: CorrelationResult, reference: Alignment):
        margin = len(result.query.sequence) * self.resolution / 2
        peaks = Peaks(result)
        return reference.refStartPosition - margin <= peaks.max * self.resolution <= reference.refEndPosition + margin
