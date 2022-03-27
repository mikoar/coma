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
                 perfectMatchScore: int,
                 scoreMultiplier: int,
                 unmatchedPenalty: int,
                 minScore: float,
                 breakSegmentThreshold: float):
        self.perfectMatchScore = perfectMatchScore
        self.scoreMultiplier = scoreMultiplier
        self.unmatchedPenalty = unmatchedPenalty
        self.minScore = minScore
        self.breakSegmentThreshold = breakSegmentThreshold

    def getSegments(self, positions: List[ScoredAlignmentPosition]) -> List[AlignmentSegment]:
        return [AlignmentSegment(positions, 0, len(positions), 0)]


class AlignmentSegments:
    @staticmethod
    def filterSegments(positions,
                       perfectMatchScore: int,
                       scoreMultiplier: int,
                       unmatchedPenalty: int,
                       minScore: float,
                       breakSegmentThreshold: float):

        scores = AlignmentSegments.getScoredPositions(positions, perfectMatchScore, scoreMultiplier, unmatchedPenalty)
        segments = AlignmentSegments.getSegments(scores, minScore, breakSegmentThreshold)
        # filteredPositions = [position for segment in segments for position in segment.positions]
        # totalScore = sum(segment.segmentScore for segment in segments)
        # return dataclasses.replace(alignment, positions=filteredPositions, confidence=totalScore)
        return segments

    @staticmethod
    def getScoredPositions(positions, perfectMatchScore: int, scoreMultiplier: int, unmatchedPenalty: int):
        return [p.getScoredPosition(perfectMatchScore, scoreMultiplier, unmatchedPenalty) for p in
                positions]

    @staticmethod
    def getSegments(positions: List[ScoredAlignmentPosition], minScore: float, breakSegmentThreshold: float) \
            -> List[AlignmentSegment]:
        if minScore <= 0:
            raise ValueError("minScore has to be bigger than 0")
        start = end = 0
        currentScore = 0
        resultSegments = []
        segmentWithMaxScore = AlignmentSegment.create(positions, 0, 0, 0)

        alignmentEnd = len(positions) - 1
        while end <= alignmentEnd:
            currentScore += positions[end].score
            if currentScore > max(0., segmentWithMaxScore.segmentScore - breakSegmentThreshold):
                end += 1
                if currentScore > segmentWithMaxScore.segmentScore:
                    segmentWithMaxScore = AlignmentSegment.create(positions, start, end, currentScore)
            else:
                if segmentWithMaxScore.segmentScore >= minScore:
                    resultSegments.append(segmentWithMaxScore)
                    segmentWithMaxScore = AlignmentSegment.create(positions, 0, 0, 0)
                start = end = end + 1
                currentScore = 0
        if segmentWithMaxScore.segmentScore >= minScore:
            resultSegments.append(segmentWithMaxScore)

        return resultSegments
