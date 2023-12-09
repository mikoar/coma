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
        self.cursor = 0
        self.currentSegmentScore = 0
        self.resultSegments = []
        self.segmentWithMaxScore = EmptyAlignmentSegment(peak, positions)

    def getSegments(self) -> List[AlignmentSegment]:
        alignmentEnd = len(self.positions) - 1
        while self.cursor <= alignmentEnd:
            self.currentSegmentScore += self.positions[self.cursor].score
            if self.__currentSegmentScoreFellBelowBreakSegmentThreshold():
                self.__breakSegment()
            else:
                self.cursor += 1
                self.__extendCurrentSegmentToCursorIfItImprovesTheScore()

        self.__addCurrentSegmentToResultIfScoreIsEnough()

        return self.resultSegments or [EmptyAlignmentSegment(self.peak, self.positions)]

    def __currentSegmentScoreFellBelowBreakSegmentThreshold(self):
        return self.currentSegmentScore <= max(0., self.segmentWithMaxScore.segmentScore - self.breakSegmentThreshold)

    def __breakSegment(self):
        self.__addCurrentSegmentToResultIfScoreIsEnough()
        self.currentSegmentStart = self.cursor = self.cursor + 1
        self.currentSegmentScore = 0

    def __addCurrentSegmentToResultIfScoreIsEnough(self):
        if self.segmentWithMaxScore.segmentScore >= self.minScore:
            self.resultSegments.append(self.segmentWithMaxScore)
            self.segmentWithMaxScore = EmptyAlignmentSegment(self.peak, self.positions)

    def __extendCurrentSegmentToCursorIfItImprovesTheScore(self):
        if self.currentSegmentScore > self.segmentWithMaxScore.segmentScore:
            self.segmentWithMaxScore = AlignmentSegment.create(self.positions[self.currentSegmentStart:self.cursor], self.currentSegmentScore,
                                                               self.peak, self.positions)
