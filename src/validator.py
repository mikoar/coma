from __future__ import annotations

from .alignment import Alignment
from .peak import Peak


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peak: Peak | None, reference: Alignment):
        if not peak:
            return False

        margin = reference.queryLength
        return reference.expectedPeakPosition - margin <= peak.positionInReference <= reference.expectedPeakPosition + margin
