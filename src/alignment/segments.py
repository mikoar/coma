from __future__ import annotations

from typing import List

from src.alignment.alignment_position import ScoredAlignmentPosition


class AlignmentSegment:
    def __init__(self, positions: List[ScoredAlignmentPosition], start: int, end: int, segmentScore: float):
        self.positions = positions
        self.start = start
        self.end = end
        self.segmentScore = segmentScore

    @staticmethod
    def create(allPositions: List[ScoredAlignmentPosition], start: int, end: int, segmentScore: float):
        positions = allPositions[start:end]
        return AlignmentSegment(positions, start, end, segmentScore)


class AlignmentSegmentsFactory:
    def __init__(self,
                 minScore: float,
                 breakSegmentThreshold: float):
        if minScore <= 0:
            raise ValueError("minScore has to be bigger than 0")
        self.minScore = minScore
        self.breakSegmentThreshold = breakSegmentThreshold

    def getSegments(self, positions: List[ScoredAlignmentPosition]) -> List[AlignmentSegment]:
        start = end = 0
        currentScore = 0
        resultSegments = []
        segmentWithMaxScore = AlignmentSegment.create(positions, 0, 0, 0)

        alignmentEnd = len(positions) - 1
        while end <= alignmentEnd:
            currentScore += positions[end].score
            if currentScore > max(0., segmentWithMaxScore.segmentScore - self.breakSegmentThreshold):
                end += 1
                if currentScore > segmentWithMaxScore.segmentScore:
                    segmentWithMaxScore = AlignmentSegment.create(positions, start, end, currentScore)
            else:
                if segmentWithMaxScore.segmentScore >= self.minScore:
                    resultSegments.append(segmentWithMaxScore)
                    segmentWithMaxScore = AlignmentSegment.create(positions, 0, 0, 0)
                start = end = end + 1
                currentScore = 0
        if segmentWithMaxScore.segmentScore >= self.minScore:
            resultSegments.append(segmentWithMaxScore)

        return resultSegments
