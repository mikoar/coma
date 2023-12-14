from typing import List

from src.alignment.alignment_position import ScoredAlignmentPosition
from src.alignment.segments import AlignmentSegment, EmptyAlignmentSegment
from src.correlation.peak import Peak


class AlignmentSegmentsFactory:
    def __init__(self,
                 minScore: float,
                 breakSegmentThreshold: float):
        if minScore <= 0:
            raise ValueError("minScore has to be bigger than 0")
        self.minScore = minScore
        self.breakSegmentThreshold = breakSegmentThreshold

    def getSegments(self, positions: List[ScoredAlignmentPosition], peak: Peak) -> List[AlignmentSegment]:
        return _AlignmentSegmentBuilder(
            self.minScore,
            self.breakSegmentThreshold,
            positions,
            peak).getSegments()


class _AlignmentSegmentBuilder:
    def __init__(self,
                 minScore: float,
                 breakSegmentThreshold: float,
                 positions: List[ScoredAlignmentPosition],
                 peak: Peak):
        self.minScore = minScore
        self.breakSegmentThreshold = breakSegmentThreshold
        self.positions = positions
        self.peak = peak
        self.currentSegmentStart = 0
        self.extendedSegmentEndPosition = 0
        self.extendedSegmentScore = 0
        self.currentSegment = EmptyAlignmentSegment(peak, positions)
        self.resultSegments = []

    def getSegments(self) -> List[AlignmentSegment]:
        alignmentEnd = len(self.positions) - 1
        while self.extendedSegmentEndPosition <= alignmentEnd:
            self.extendedSegmentScore += self.positions[self.extendedSegmentEndPosition].score
            if self.__extendedSegmentScoreFellBelowBreakSegmentThreshold():
                self.__breakSegment()
            else:
                self.extendedSegmentEndPosition += 1
                self.__acceptExtendedSegmentIfScoreIsImproved()

        self.__addCurrentSegmentToResultIfScoreIsEnough()

        return self.resultSegments or [EmptyAlignmentSegment(self.peak, self.positions)]

    def __extendedSegmentScoreFellBelowBreakSegmentThreshold(self):
        return self.extendedSegmentScore <= max(0., self.currentSegment.segmentScore - self.breakSegmentThreshold)

    def __breakSegment(self):
        self.__addCurrentSegmentToResultIfScoreIsEnough()
        self.currentSegmentStart = self.extendedSegmentEndPosition = self.extendedSegmentEndPosition + 1
        self.extendedSegmentScore = 0

    def __addCurrentSegmentToResultIfScoreIsEnough(self):
        if self.currentSegment.segmentScore >= self.minScore:
            self.resultSegments.append(self.currentSegment)
            self.currentSegment = EmptyAlignmentSegment(self.peak, self.positions)

    def __acceptExtendedSegmentIfScoreIsImproved(self):
        if self.extendedSegmentScore > self.currentSegment.segmentScore:
            self.currentSegment = AlignmentSegment.create(
                self.positions[self.currentSegmentStart:self.extendedSegmentEndPosition],
                self.peak,
                self.positions)
