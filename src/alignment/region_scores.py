from __future__ import annotations

from typing import List, NamedTuple


class AlignmentSegment(NamedTuple):
    start: int
    end: int
    score: float


class RegionScores:
    def __init__(self, scores: List[float]):
        self.scores = scores

    def getSegments(self, minScore, breakSegmentThreshold) -> List[AlignmentSegment]:
        start = end = 0
        currentScore = 0
        resultSegments = []
        segmentWithMaxScore = AlignmentSegment(0, 0, 0)

        alignmentEnd = len(self.scores) - 1
        while end <= alignmentEnd:
            currentScore += self.scores[end]
            if currentScore > max(0, segmentWithMaxScore.score - breakSegmentThreshold):
                if currentScore > segmentWithMaxScore.score:
                    segmentWithMaxScore = AlignmentSegment(start, end, currentScore)
                end += 1
            else:
                if segmentWithMaxScore.score >= minScore:
                    resultSegments.append(segmentWithMaxScore)
                    segmentWithMaxScore = AlignmentSegment(0, 0, 0)
                start = end = end + 1
                currentScore = 0
        if segmentWithMaxScore.score >= minScore:
            resultSegments.append(segmentWithMaxScore)

        return resultSegments
