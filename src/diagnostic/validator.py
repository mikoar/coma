from __future__ import annotations

from src.correlation.bionano_alignment import BionanoAlignment
from src.correlation.peak import Peak


class Validator:
    def __init__(self, resolution: int, tolerance: int = 1024) -> None:
        self.resolution = resolution
        self.tolerance = tolerance

    def validate(self, peak: Peak | None, reference: BionanoAlignment):
        if not peak:
            return False

        return self.__peakWithinAlignmentSizeUncertainty(peak, reference)

    def __peakWithinAlignmentSizeUncertainty(self, peak: Peak, reference: BionanoAlignment):
        margin = abs(reference.queryReferenceAlignmentLengthDifference) / 2 + self.tolerance
        return reference.referenceStartPosition - margin <= peak.position <= reference.referenceStartPosition + margin
