from __future__ import annotations

from src.correlation.bionano_alignment import BionanoAlignment
from src.correlation.peak import Peak


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peak: Peak | None, reference: BionanoAlignment):
        if not peak:
            return False

        return self.__peakWithinAlignmentSizeUncertaintyFromCenterWithFixedMargin(peak, reference)

    def __peakWithinAlignmentSizeFromCenter(self, peak: Peak, reference: BionanoAlignment):
        margin = max(reference.queryAlignmentLength(), reference.referenceAlignmentLength()) / 2 + self.resolution
        return reference.expectedPeakPosition - margin <= peak.position <= reference.expectedPeakPosition + margin

    def __peakWithinAlignmentSizeUncertaintyFromCenter(self, peak: Peak, reference: BionanoAlignment):
        margin = abs(reference.queryReferenceAlignmentLengthDifference) / 2 + self.resolution
        return reference.expectedPeakPosition - margin <= peak.position <= reference.expectedPeakPosition + margin

    def __peakWithinAlignmentSizeUncertaintyFromCenterWithFixedMargin(self, peak: Peak, reference: BionanoAlignment):
        margin = abs(reference.queryReferenceAlignmentLengthDifference) / 2 + 1024
        return reference.expectedPeakPosition - margin <= peak.position <= reference.expectedPeakPosition + margin

    def __peakAnywhereInMolecule(self, peak: Peak, reference: BionanoAlignment):
        return reference.expectedQueryMoleculeStart - self.resolution <= peak.position <= reference.expectedQueryMoleculeEnd + self.resolution
