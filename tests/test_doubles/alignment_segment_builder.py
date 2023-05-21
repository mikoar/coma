from typing import List

from src.alignment.alignment_position import ScoredAlignmentPosition
from src.alignment.segments import AlignmentSegment
from src.correlation.peak import Peak


class AlignmentSegmentBuilder:
    def __init__(self):
        self.positions: List[ScoredAlignmentPosition] = []
        self.segmentScore = 1200.
        self.peak: Peak = Peak(0, 0.)

    def withPosition(self, position: ScoredAlignmentPosition):
        self.positions.append(position)
        return self

    def withScore(self, segmentScore: float):
        self.segmentScore = segmentScore
        return self

    def build(self):
        return AlignmentSegment(self.positions, self.segmentScore, self.peak)
