from __future__ import annotations

import itertools
from typing import List

from src.alignment.alignment_position import ScoredAlignmentPosition, ScoredAlignedPair, AlignedPair


class AlignmentSegment:
    def __init__(self, positions: List[ScoredAlignmentPosition], segmentScore: float):
        self.positions = positions
        self.segmentScore = segmentScore
        self.alignedPositions = [p for p in positions if isinstance(p, ScoredAlignedPair)]
        # self.referencePositions = [p.reference for p in self.alignedPositions]
        # self.queryPositions = [p.query for p in self.alignedPositions]

    def sliceByReference(self, startSiteId: int, endSiteId) -> AlignmentSegment:
        slicedAtStart = itertools.dropwhile(
            lambda p: not isinstance(p, AlignedPair) or p.reference.siteId < startSiteId,
            self.positions)
        positions = list(itertools.takewhile(lambda p: not isinstance(p, AlignedPair) or p.reference.siteId < endSiteId,
                                             slicedAtStart))
        self.trimNotAlignedPositions(positions)
        return AlignmentSegment(positions, sum(p.score for p in positions))

    @staticmethod
    def trimNotAlignedPositions(positions):
        while not isinstance(positions[-1], AlignedPair):
            positions.pop()

    def getConflict(self, other: AlignmentSegment):
        otherReferenceStart = other.alignedPositions[0].reference
        otherReferenceEnd = other.alignedPositions[-1].reference
        conflictingReferencePositions = [p for p in self.alignedPositions if
                                         otherReferenceStart.siteId <= p.reference.siteId <= otherReferenceEnd.siteId]


class __EmptyAlignmentSegment(AlignmentSegment):
    def __init__(self):
        super().__init__([], 0.)


AlignmentSegment.empty = __EmptyAlignmentSegment()


class _Conflict:
    def __init__(self, segment1: AlignmentSegment, conflictingPositions1: List[ScoredAlignmentPosition],
                 segment2: AlignmentSegment, conflictingPositions2: List[ScoredAlignmentPosition]):
        self.segment1 = segment1
        self.conflictingPositions1 = conflictingPositions1
        self.segment2 = segment2
        self.conflictingPositions2 = conflictingPositions2

    # @staticmethod
    # def referenceConflict(segment1: AlignmentSegment, segment2: AlignmentSegment,
    #                       conflictStartSiteId: int, conflictEndSiteId: int):
    #     return _Conflict(segment1, [p for p in segment1.positions if p.])


class AlignmentSegmentsWithoutConflicts:
    @staticmethod
    def create(segments: List[AlignmentSegment]):
        # refDupes = itertools.groupby(positionsWithSegments, lambda x: x.position.)
        return AlignmentSegmentsWithoutConflicts(segments)

    def __init__(self, segments: List[AlignmentSegment]):
        self.segments = segments


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
        segmentWithMaxScore = AlignmentSegment.empty

        alignmentEnd = len(positions) - 1
        while end <= alignmentEnd:
            currentScore += positions[end].score
            if currentScore > max(0., segmentWithMaxScore.segmentScore - self.breakSegmentThreshold):
                end += 1
                if currentScore > segmentWithMaxScore.segmentScore:
                    segmentWithMaxScore = AlignmentSegment(positions[start:end], currentScore)
            else:
                if segmentWithMaxScore.segmentScore >= self.minScore:
                    resultSegments.append(segmentWithMaxScore)
                    segmentWithMaxScore = AlignmentSegment.empty
                start = end = end + 1
                currentScore = 0
        if segmentWithMaxScore.segmentScore >= self.minScore:
            resultSegments.append(segmentWithMaxScore)

        return resultSegments
