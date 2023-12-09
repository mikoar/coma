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
               segmentScore: float,
               peak: Peak,
               allPeakPositions: List[ScoredAlignmentPosition]):
        return AlignmentSegment(positions, segmentScore, peak, allPeakPositions) \
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
            return _SegmentPairWithConflict.create(self, other, other.startPosition, self.endPosition)
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
        return AlignmentSegment.create(positions, sum(p.score for p in positions), self.peak, self.allPeakPositions)

    def getReferenceLabels(self) -> _SegementMoleculeCharacteristics:
        """Function used to get all of the reference labels,
        scores and their indexes present in a segment

        :return: Characteristics of Reference in a segment of alignment
        :rtype: _SegementMoleculeCharacteristics
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
        return _SegementMoleculeCharacteristics(referencePositions, referenceScores, referenceIndexes)

    def getQueryLabels(self) -> _SegementMoleculeCharacteristics:
        """Function used to get all of the query labels,
        scores and their indexes present in a segment

        :return: Characteristics of Query in a segment of alignment
        :rtype: _SegementMoleculeCharacteristics
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
        return _SegementMoleculeCharacteristics(queryPositions, queryScores, queryIndexes)

    @staticmethod
    def __trimNotAlignedPositionsFromEnd(positions, end=None):
        if positions:
            while not isinstance(positions[-1], AlignedPair) and (
                    not end or not positions[-1].lessOrEqualOnAnySequence(end)):
                positions.pop()

    def __eq__(self, other):
        return isinstance(other, AlignmentSegment) \
               and other.segmentScore == self.segmentScore \
               and other.positions == self.positions \
               and other.peak == self.peak

    def __sub__(self, other: AlignmentSegment):
        positions = [p for p in self.positions if p not in other.positions]
        return AlignmentSegment.create(positions, sum(p.score for p in positions), self.peak, self.allPeakPositions)

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
    def __init__(self, segment1: AlignmentSegment, segment2: AlignmentSegment):
        self.segment1 = segment1
        self.segment2 = segment2

    @abstractmethod
    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        pass


class _SegmentPairWithNoConflict(_SegmentPair):
    def __init__(self, segment1: AlignmentSegment, segment2: AlignmentSegment):
        super().__init__(segment1, segment2)

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        return self.segment1, self.segment2


class _SegmentPairWithConflict(_SegmentPair):
    def __init__(self, segment1: AlignmentSegment, conflictingSubsegment1: AlignmentSegment,
                 segment2: AlignmentSegment, conflictingSubsegment2: AlignmentSegment):
        super().__init__(segment1, segment2)
        self.conflictingSubsegment1 = conflictingSubsegment1
        self.conflictingSubsegment2 = conflictingSubsegment2

    @staticmethod
    def create(segment1: AlignmentSegment, segment2: AlignmentSegment,
               conflictStart: AlignedPair, conflictEnd: AlignedPair):
        conflictingSubsegment1 = segment1.slice(conflictStart, conflictEnd)
        conflictingSubsegment2 = segment2.slice(conflictStart, conflictEnd)
        return _SegmentPairWithConflict(segment1, conflictingSubsegment1, segment2, conflictingSubsegment2)

    def __findOptimaPlace(self, conf1, conf2) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if len(conf1.positions) == len(conf2.positions):
            points_1_np = np.cumsum([0] + conf1.scores)
            points_2_np = np.cumsum([0] + conf2.scores[::-1])[::-1]
            scores_sum = np.add(points_1_np, points_2_np)
            max_index = np.argmax(scores_sum)
            if max_index == 0:
                return self.segment1 - self.conflictingSubsegment1, self.segment2
            elif max_index == len(conf1.positions):
                return self.segment1, self.segment2 - self.conflictingSubsegment2
            else:
                new_seg1 = self.segment1 - AlignmentSegment.create(
                    self.conflictingSubsegment1.positions[conf1.indexes[max_index]:],
                    sum(p.score for p in self.conflictingSubsegment1.positions[conf1.indexes[max_index]:]),
                    self.segment1.peak,
                    self.conflictingSubsegment1.positions[conf1.indexes[max_index]:])
                new_seg2 = self.segment2 - AlignmentSegment.create(
                    self.conflictingSubsegment2.positions[:conf2.indexes[max_index]],
                    sum(p.score for p in self.conflictingSubsegment2.positions[:conf2.indexes[max_index]]),
                    self.segment2.peak,
                    self.conflictingSubsegment2.positions[:conf2.indexes[max_index]])
            return new_seg1, new_seg2
        else:
            self.resolveByTrimming()

    def resolveByTrimming(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        """Function used to resolve conflicts
        with uneven number of conflicting labels

        :return: Two Segments without conflicts
        :rtype: Tuple[AlignmentSegment, AlignmentSegment]
        """
        if self.conflictingSubsegment1.segmentScore > self.conflictingSubsegment2.segmentScore:
            return self.segment1, self.segment2 - self.conflictingSubsegment2
        else:
            return self.segment1 - self.conflictingSubsegment1, self.segment2

    def resolveConflict(self) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if self.conflictingSubsegment1.peak.position > self.conflictingSubsegment2.peak.position:
            return self.__findOptimaPlace(self.conflictingSubsegment1.getReferenceLabels(),
                                          self.conflictingSubsegment2.getReferenceLabels())
        else:
            return self.__findOptimaPlace(self.conflictingSubsegment1.getQueryLabels(),
                                          self.conflictingSubsegment2.getQueryLabels())


class _SegementMoleculeCharacteristics():
    def __init__(self, positions: List, scores: List, indexes: List):
        self.positions = positions
        self.scores = scores
        self.indexes = indexes

    @staticmethod
    def create(positions: List, scores: List, indexes: List):
        return _SegementMoleculeCharacteristics(positions, scores, indexes)
