from src.alignment.alignment_position import AlignedPair
from src.correlation.optical_map import PositionWithSiteId


class AlignedPairBuilder:
    def __init__(self):
        self.referenceSiteId = 1
        self.referencePosition = 0
        self.querySiteId = 1
        self.queryPosition = 0
        self.queryShift = 0

    def withReferenceSiteId(self, referenceSiteId):
        self.referenceSiteId = referenceSiteId
        return self

    def withReferencePosition(self, referencePosition):
        self.referencePosition = referencePosition
        return self

    def withQuerySiteId(self, querySiteId):
        self.querySiteId = querySiteId
        return self

    def withQueryPosition(self, queryPosition):
        self.queryPosition = queryPosition
        return self

    def withQueryShift(self, queryShift):
        self.queryShift = queryShift
        return self

    def build(self):
        return AlignedPair(
            PositionWithSiteId(self.referenceSiteId, self.referencePosition),
            PositionWithSiteId(self.querySiteId, self.queryPosition),
            self.queryShift)
