from __future__ import annotations

import itertools
from abc import ABC, abstractmethod
from typing import List, Tuple

import numpy as np

from src.alignment.alignment_position import ScoredAlignmentPosition, ScoredAlignedPair, AlignedPair, \
    NotAlignedQueryPosition, ScoredNotAlignedPosition, NotAlignedReferencePosition
from src.correlation.peak import Peak


class AlignmentSegment:
    @staticmethod
    def create(positions: List[ScoredAlignmentPosition],
               peak: Peak,
               allPeakPositions: List[ScoredAlignmentPosition]):
        return AlignmentSegment(positions, sum(p.score for p in positions), peak, allPeakPositions) \
            if positions \
            else EmptyAlignmentSegment(peak, allPeakPositions)

    def __init__(
            self,
            positions: List[ScoredAlignmentPosition],
            segmentScore: float,
            peak: Peak,
            allPeakPositions: List[ScoredAlignmentPosition]):
        self.positions = positions
        self.segmentScore = segmentScore
        self.alignedPositions = [p for p in positions if isinstance(p, ScoredAlignedPair)]
        self.peak = peak
        self.allPeakPositions = allPeakPositions or []

    @property
    def empty(self):
        return len(self.positions) == 0

    @property
    def startPosition(self):
        return self.alignedPositions[0]

    @property
    def endPosition(self):
        return self.alignedPositions[-1]

    @property
    def reverse(self):
        return self.startPosition.query.siteId > self.endPosition.query.siteId

    def checkForConflicts(self, other: AlignmentSegment):
        if self.endOverlapsWithStartOf(other):
            return _SegmentPairWithConflict.create(self, other)
        return _SegmentPairWithNoConflict(self, other)

    def endOverlapsWithStartOf(self, other: AlignmentSegment):
        return other.startPosition.lessOrEqualOnAnySequence(self.startPosition) or \
               other.startPosition.lessOrEqualOnAnySequence(self.endPosition) or \
               self.endPosition.lessOrEqualOnAnySequence(other.endPosition)

    def slice(self, start: AlignedPair, end: AlignedPair) -> AlignmentSegment:
        slicedAtStart = itertools.dropwhile(
            lambda p: p.lessOnBothSequences(start), self.positions)
        positions = list(
            itertools.takewhile(
                lambda p: not isinstance(p, AlignedPair) or p.lessOrEqualOnAnySequence(end),
                slicedAtStart))
        self.__trimNotAlignedPositionsFromEnd(positions, end)
        return AlignmentSegment.create(positions, self.peak, self.allPeakPositions)

    def getReferenceLabels(self) -> _ConflictingSegmentCharacteristics:
        """Function used to get all of the reference labels,
        scores and their indexes present in a segment

        :return: Characteristics of Reference in a segment of alignment
        :rtype: _ConflictingSegmentCharacteristics
        """

        referencePositions, referenceScores, referenceIndexes = [], [], []
        sumScore = 0
        for index, position in enumerate(self.positions):
            if isinstance(position, ScoredAlignedPair):
                referenceIndexes.append(index)
                referencePositions.append(position.reference)
                referenceScores.append(position.score + sumScore)
                sumScore = 0
            elif isinstance(position, ScoredNotAlignedPosition) \
                    and isinstance(position.position, NotAlignedReferencePosition):
                referenceIndexes.append(index)
                referencePositions.append(position.position.reference)
                referenceScores.append(position.score + sumScore)
                sumScore = 0
            else:
                sumScore += position.score
        return _ConflictingSegmentCharacteristics(referencePositions, referenceScores, referenceIndexes)

    def getQueryLabels(self) -> _ConflictingSegmentCharacteristics:
        """Function used to get all of the query labels,
        scores and their indexes present in a segment

        :return: Characteristics of Query in a segment of alignment
        :rtype: _ConflictingSegmentCharacteristics
        """
        queryPositions, queryScores, queryIndexes = [], [], []
        sumScore = 0
        for index, position in enumerate(self.positions):
            if isinstance(position, ScoredAlignedPair):
                queryPositions.append(position.query)
                queryScores.append(position.score + sumScore)
                queryIndexes.append(index)
                sumScore = 0
            elif isinstance(position, ScoredNotAlignedPosition) \
                    and isinstance(position.position, NotAlignedQueryPosition):
                queryPositions.append(position.position.query)
                queryScores.append(position.score + sumScore)
                queryIndexes.append(index)
                sumScore = 0
            else:
                sumScore += position.score
        return _ConflictingSegmentCharacteristics(queryPositions, queryScores, queryIndexes)

    @staticmethod
    def __trimNotAlignedPositionsFromEnd(positions, end=None):
        while positions and (not isinstance(positions[-1], AlignedPair)) and (
                    not end or not positions[-1].lessOrEqualOnAnySequence(end)):
            positions.pop()
        # moze skonczyc sie pustym fragmentem - do sprawdzenia
        # if positions:
        #     while not isinstance(positions[-1], AlignedPair) and (
        #             not end or not positions[-1].lessOrEqualOnAnySequence(end)):
        #         positions.pop()

    def __eq__(self, other):
        return isinstance(other, AlignmentSegment) \
               and other.segmentScore == self.segmentScore \
               and other.positions == self.positions \
               and other.peak == self.peak

    def __sub__(self, other: AlignmentSegment | List[ScoredAlignmentPosition]):
        if isinstance(other, AlignmentSegment):
            otherPositions = other.positions
        else:
            otherPositions = other
        positions = [p for p in self.positions if p not in otherPositions]
        return AlignmentSegment.create(positions, self.peak, self.allPeakPositions)

    def __repr__(self):
        return f"score: {self.segmentScore}, positions: {self.positions}"


class EmptyAlignmentSegment(AlignmentSegment):
    def __init__(self,
                 peak: Peak = None,
                 allPeakPositions: List[ScoredAlignmentPosition] = None):
        super().__init__([], 0., peak or Peak.null, allPeakPositions or [])

    @property
    def startPosition(self):
        return AlignedPair.null

    @property
    def endPosition(self):
        return AlignedPair.null

    def checkForConflicts(self, other: AlignmentSegment):
        return _SegmentPairWithNoConflict(self, other)

    def endOverlapsWithStartOf(self, other: AlignmentSegment):
        return False


class _SegmentPair(ABC):
    def __init__(self, leftSegment: AlignmentSegment, rightSegment: AlignmentSegment):
        self.leftSegment = leftSegment
        self.rightSegment = rightSegment

    @abstractmethod
    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        pass


class _SegmentPairWithNoConflict(_SegmentPair):
    def __init__(self, leftSegment: AlignmentSegment, rightSegment: AlignmentSegment):
        super().__init__(leftSegment, rightSegment)

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        return self.leftSegment, self.rightSegment


class _SegmentPairWithConflict(_SegmentPair):
    def __init__(self, leftSegment: AlignmentSegment, leftConflictingSubsegment: AlignmentSegment,
                 rightSegment: AlignmentSegment, rightConflictingSubsegment: AlignmentSegment):
        super().__init__(leftSegment, rightSegment)
        self.leftConflictingSubsegment = leftConflictingSubsegment
        self.rightConflictingSubsegment = rightConflictingSubsegment

    @staticmethod
    def create(segment1: AlignmentSegment, segment2: AlignmentSegment):
        conflictStart = segment2.startPosition
        conflictEnd = segment1.endPosition
        leftConflictingSubsegment = segment1.slice(conflictStart, conflictEnd)
        rightConflictingSubsegment = segment2.slice(conflictStart, conflictEnd)
        return _SegmentPairWithConflict(segment1, leftConflictingSubsegment, segment2, rightConflictingSubsegment)

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if self.leftConflictingSubsegment.peak.position > self.rightConflictingSubsegment.peak.position:
            return self.__trimSegmentsAtOptimalPosition(
                self.leftConflictingSubsegment.getReferenceLabels(),
                self.rightConflictingSubsegment.getReferenceLabels())
        else:
            return self.__trimSegmentsAtOptimalPosition(
                self.leftConflictingSubsegment.getQueryLabels(),
                self.rightConflictingSubsegment.getQueryLabels())

    def __trimSegmentsAtOptimalPosition(
            self,
            leftSubsegmentCharacteristics: _ConflictingSegmentCharacteristics,
            rightSubsegmentCharacteristics: _ConflictingSegmentCharacteristics):
        if len(leftSubsegmentCharacteristics.positions) == len(rightSubsegmentCharacteristics.positions):
            optimalMergeIndex = self.__getOptimalMergeIndex(leftSubsegmentCharacteristics,
                                                            rightSubsegmentCharacteristics)
            if optimalMergeIndex == 0:
                return self.leftSegment - self.leftConflictingSubsegment, self.rightSegment
            elif optimalMergeIndex == len(leftSubsegmentCharacteristics.positions):
                return self.leftSegment, self.rightSegment - self.rightConflictingSubsegment
            else:
                leftTrimIndex = leftSubsegmentCharacteristics.indexes[optimalMergeIndex]
                leftSegmentPositionsToRemove = self.leftConflictingSubsegment.positions[leftTrimIndex:]
                newLeftSegment = self.leftSegment - leftSegmentPositionsToRemove

                rightTrimIndex = rightSubsegmentCharacteristics.indexes[optimalMergeIndex]
                rightSegmentPositionsToRemove = self.rightConflictingSubsegment.positions[:rightTrimIndex]
                newRightSegment = self.rightSegment - rightSegmentPositionsToRemove
            return newLeftSegment, newRightSegment
        else:
            return self.__removeWholeConflictingSubsegmentWithWorseScore()

    @staticmethod
    def __getOptimalMergeIndex(leftSubsegmentCharacteristics, rightSubsegmentCharacteristics):
        leftSubsegmentCumulatedScores = np.cumsum([0] + leftSubsegmentCharacteristics.scores)
        rightSubsegmentCumulatedScores = np.cumsum([0] + rightSubsegmentCharacteristics.scores[::-1])[::-1]
        totalCumulatedScores = np.add(leftSubsegmentCumulatedScores, rightSubsegmentCumulatedScores)
        optimalMergeIndex = np.argmax(totalCumulatedScores)
        return optimalMergeIndex

    def __removeWholeConflictingSubsegmentWithWorseScore(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if self.leftConflictingSubsegment.segmentScore > self.rightConflictingSubsegment.segmentScore:
            return self.leftSegment, self.rightSegment - self.rightConflictingSubsegment
        else:
            return self.leftSegment - self.leftConflictingSubsegment, self.rightSegment


class _ConflictingSegmentCharacteristics:
    def __init__(self, positions: List, scores: List, indexes: List):
        self.positions = positions
        self.scores = scores
        self.indexes = indexes
