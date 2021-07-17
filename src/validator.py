from cmap_reader import Alignment
from optical_map import Peaks


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peaks: Peaks, reference: Alignment):
        start = reference.refStartPosition - reference.queryStartPosition
        end = reference.refEndPosition + peaks.queryLength - reference.queryEndPosition
        return start <= peaks.max * self.resolution <= end
