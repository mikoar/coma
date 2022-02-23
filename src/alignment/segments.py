from __future__ import annotations

import dataclasses
from typing import List

from src.alignment.alignment_position import ScoredAlignmentPosition
from src.alignment.alignment_results import AlignmentResultRow


class AlignmentSegment:
    def __init__(self, allPositions: List[ScoredAlignmentPosition], start: int, end: int, segmentScore: float):
        self.positions = allPositions[start:end]
        self.start = start
        self.end = end
        self.segmentScore = segmentScore


class AlignmentSegments:
    @staticmethod
    def filterSegments(alignment: AlignmentResultRow,
                       perfectMatchScore: int,
                       scoreMultiplier: int,
                       unmatchedPenalty: int,
                       minScore: float,
                       breakSegmentThreshold: float):

        scores = alignment.getScoredPositions(perfectMatchScore, scoreMultiplier, unmatchedPenalty)
        segments = AlignmentSegments.getSegments(scores, minScore, breakSegmentThreshold)
        filteredPositions = [position for segment in segments for position in segment.positions]
        totalScore = sum(segment.segmentScore for segment in segments)
        return dataclasses.replace(alignment, positions=filteredPositions, confidence=totalScore)

    @staticmethod
    def getSegments(positions: List[ScoredAlignmentPosition], minScore: float, breakSegmentThreshold: float) \
            -> List[AlignmentSegment]:
        if minScore <= 0:
            raise ValueError("minScore has to be bigger than 0")
        start = end = 0
        currentScore = 0
        resultSegments = []
        segmentWithMaxScore = AlignmentSegment(positions, 0, 0, 0)

        alignmentEnd = len(positions) - 1
        while end <= alignmentEnd:
            currentScore += positions[end].score
            if currentScore > max(0., segmentWithMaxScore.segmentScore - breakSegmentThreshold):
                end += 1
                if currentScore > segmentWithMaxScore.segmentScore:
                    segmentWithMaxScore = AlignmentSegment(positions, start, end, currentScore)
            else:
                if segmentWithMaxScore.segmentScore >= minScore:
                    resultSegments.append(segmentWithMaxScore)
                    segmentWithMaxScore = AlignmentSegment(positions, 0, 0, 0)
                start = end = end + 1
                currentScore = 0
        if segmentWithMaxScore.segmentScore >= minScore:
            resultSegments.append(segmentWithMaxScore)

        return resultSegments
