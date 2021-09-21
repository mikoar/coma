from typing import List, NamedTuple


class Segment(NamedTuple):
    start: int
    end: int
    segmentScore: float


class SegmentScorer:
    def getSegmentWithMaxScore(self, alignmentScores: List[float]):
        start = end = 0
        currentScore = 0.
        maxScore = Segment(0, 0, 0)

        alignmentEnd = len(alignmentScores) - 1
        while end <= alignmentEnd:
            if currentScore + alignmentScores[end] > 0:
                currentScore += alignmentScores[end]
                if currentScore > maxScore.segmentScore:
                    maxScore = Segment(start, end, currentScore)
                end += 1
            else:
                start = end = end + 1
                currentScore = 0

        return maxScore
