from typing import List, NamedTuple


class AlignmentSegment(NamedTuple):
    start: int
    end: int
    score: float


class AlignmentScorer:
    def getSegmentWithMaxScore(self, alignmentScores: List[float]):
        start = end = 0
        currentScore = 0
        segmentWithMaxScore = AlignmentSegment(0, 0, 0)

        alignmentEnd = len(alignmentScores) - 1
        while end <= alignmentEnd:
            if currentScore + alignmentScores[end] > 0:
                currentScore += alignmentScores[end]
                if currentScore > segmentWithMaxScore.score:
                    segmentWithMaxScore = AlignmentSegment(start, end, currentScore)
                end += 1
            else:
                start = end = end + 1
                currentScore = 0

        return segmentWithMaxScore
