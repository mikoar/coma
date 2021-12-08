from __future__ import annotations

from typing import List, NamedTuple


class AlignmentSegment(NamedTuple):
    start: int
    end: int
    score: float


class RegionScores:
    def __init__(self, scores: List[float]):
        self.scores = scores

    def getSegmentWithMaxScore(self):
        start = end = 0
        currentScore = 0
        segmentWithMaxScore = AlignmentSegment(0, 0, 0)

        alignmentEnd = len(self.scores) - 1
        while end <= alignmentEnd:
            if currentScore + self.scores[end] > 0:
                currentScore += self.scores[end]
                if currentScore > segmentWithMaxScore.score:
                    segmentWithMaxScore = AlignmentSegment(start, end, currentScore)
                end += 1
            else:
                start = end = end + 1
                currentScore = 0

        return segmentWithMaxScore

# TODO choose multiple segments
