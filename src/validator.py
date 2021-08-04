from __future__ import annotations

from .alignment import Alignment
from .peak import Peak


class Validator:
    def __init__(self, resolution: int) -> None:
        self.resolution = resolution

    def validate(self, peak: Peak | None, reference: Alignment):
        if not peak:
            return False

        return reference.expectedQueryMoleculeStart + self.resolution <= peak.positionInReference <= reference.expectedQueryMoleculeEnd + self.resolution
