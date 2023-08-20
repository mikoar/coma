import numpy as np

from src.correlation.optical_map import InitialAlignment, OpticalMap


class InitialAlignmentBuilder:
    def __init__(self):
        self.peaks = []

    def withPeak(self, peak):
        self.peaks.append(peak)
        return self

    def build(self):
        return InitialAlignment(
            np.array([]),
            OpticalMap(1, 1, []),
            OpticalMap(1, 1, []),
            self.peaks,
            False,
            0.)
