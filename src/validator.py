from __future__ import annotations

from .alignment import Alignment
from .peak import Peak


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peak: Peak | None, reference: Alignment):
        if not peak:
            return False

        return self.__x(peak, reference)

    def __x(self, peak: Peak, reference: Alignment):
        margin = max(reference.queryAlignmentLength(), reference.referenceAlignmentLength())/2 + self.resolution
        return reference.expectedPeakPosition - margin <= peak.positionInReference <= reference.expectedPeakPosition + margin

    def __peakInMoleculeCenter_evenMoreStrict(self, peak: Peak, reference: Alignment):
        margin = abs(reference.queryReferenceAlignmentLengthDifference)/2 + self.resolution
        return reference.expectedPeakPosition - margin <= peak.positionInReference <= reference.expectedPeakPosition + margin

    def __peakInMoleculeCenter(self, peak: Peak, reference: Alignment):
        margin = abs(reference.queryReferenceAlignmentLengthDifference) + self.resolution
        return reference.expectedPeakPosition - margin <= peak.positionInReference <= reference.expectedPeakPosition + margin

    def __peakAnywhereInMolecule(self, peak: Peak, reference: Alignment):
        return reference.expectedQueryMoleculeStart + self.resolution <= peak.positionInReference <= reference.expectedQueryMoleculeEnd + self.resolution
