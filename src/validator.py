from __future__ import annotations

from .alignment import Alignment
from .peak import Peak


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peak: Peak | None, reference: Alignment, restrictive=False):
        if not peak:
            return False

        return self.peakInMoleculeCenter(peak, reference) if restrictive else self.peakAnywhereInMolecule(peak, reference)

    def peakInMoleculeCenter(self, peak: Peak, reference: Alignment):
        margin = abs(reference.queryReferenceAlignmentLengthDifference) + self.resolution
        return reference.expectedPeakPosition - margin <= peak.positionInReference <= reference.expectedPeakPosition + margin

    def peakAnywhereInMolecule(self, peak: Peak, reference: Alignment):
        return reference.expectedQueryMoleculeStart + self.resolution <= peak.positionInReference <= reference.expectedQueryMoleculeEnd + self.resolution
