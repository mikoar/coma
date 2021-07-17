from cmap_reader import Alignment
from optical_map import Peaks


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peaks: Peaks, reference: Alignment):
        margin = len(peaks.correlationResult.query.sequence) * self.resolution / 2
        return reference.refStartPosition - margin <= peaks.max * self.resolution <= reference.refEndPosition + margin
