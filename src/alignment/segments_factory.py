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
        start = end = 0
        currentScore = 0
        resultSegments = []
        segmentWithMaxScore = EmptyAlignmentSegment(peak, positions)

        alignmentEnd = len(positions) - 1
        while end <= alignmentEnd:
            currentScore += positions[end].score
            if currentScore > max(0., segmentWithMaxScore.segmentScore - self.breakSegmentThreshold):
                end += 1
                if currentScore > segmentWithMaxScore.segmentScore:
                    segmentWithMaxScore = AlignmentSegment.create(positions[start:end], currentScore, peak, positions)
            else:
                if segmentWithMaxScore.segmentScore >= self.minScore:
                    resultSegments.append(segmentWithMaxScore)
                    segmentWithMaxScore = EmptyAlignmentSegment(peak, positions)
                start = end = end + 1
                currentScore = 0
        if segmentWithMaxScore.segmentScore >= self.minScore:
            resultSegments.append(segmentWithMaxScore)

        return resultSegments or [EmptyAlignmentSegment(peak, positions)]
