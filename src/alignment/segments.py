from __future__ import annotations

import itertools
import numpy as np
from abc import ABC, abstractmethod
from typing import List, Tuple

from src.alignment.alignment_position import ScoredAlignmentPosition, ScoredAlignedPair, AlignedPair, NotAlignedQueryPosition
from src.correlation.peak import Peak


class AlignmentSegment:
    empty: AlignmentSegment

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
        return AlignmentSegment(positions, sum(p.score for p in positions), self.peak, self.allPeakPositions)

    def get_reference_labels(self) -> Tuple[List, List, List]:
        """Function used to get all of the reference labels,
        scores and their indexes present in a segment

        :return: Characteristics of Reference in a segment of alignment
        :rtype: Tuple[List, List, List]
        """
        ref_pos, ref_scores, ref_index = [], [], []
        sum_score = 0
        for index, position in enumerate(self.positions):
            if isinstance(position, ScoredAlignedPair):
                ref_index.append(index)
                ref_pos.append(position.reference)
                ref_scores.append(position.score + sum_score)
                sum_score = 0
            else:
                if not isinstance(position.position, NotAlignedQueryPosition):
                    ref_index.append(index)
                    ref_pos.append(position.position.reference)
                    ref_scores.append(position.score + sum_score)
                    sum_score = 0
                else:
                    sum_score += position.score
        return (ref_pos, ref_scores, ref_index)

    def get_query_labels(self) -> Tuple(List, List, List):
        """Function used to get all of the query labels,
        scores and their indexes present in a segment

        :return: Characteristics of Query in a segment of alignment
        :rtype: Tuple[List, List, List]
        """
        quer_pos, quer_scores, quer_index = [], [], []
        sum_score = 0
        for index, position in enumerate(self.positions):
            if isinstance(position, ScoredAlignedPair):
                quer_pos.append(position.query)
                quer_scores.append(position.score + sum_score)
                quer_index.append(index)
                sum_score = 0
            else:
                if isinstance(position.position, NotAlignedQueryPosition):
                    quer_pos.append(position.position.query)
                    quer_scores.append(position.score + sum_score)
                    quer_index.append(index)
                    sum_score = 0
                else:
                    sum_score += position.score
        return (quer_pos, quer_scores, quer_index)

    @staticmethod
    def __trimNotAlignedPositionsFromEnd(positions, end=None):
        if positions:
            if not end:
                while not isinstance(positions[-1], AlignedPair):
                    positions.pop()
            else:
                while not isinstance(positions[-1], AlignedPair):
                    if not positions[-1].lessOrEqualOnAnySequence(end):
                        positions.pop()
                    else:
                        break

    def __eq__(self, other):
        return isinstance(other, AlignmentSegment) \
            and other.segmentScore == self.segmentScore \
            and other.positions == self.positions

    def __sub__(self, other: AlignmentSegment):
        positions = [p for p in self.positions if p not in other.positions]
        if not positions:
            return AlignmentSegment.empty
        return AlignmentSegment(positions, sum(p.score for p in positions), self.peak, self.allPeakPositions)

    def __repr__(self):
        return f"score: {self.segmentScore}, positions: {self.positions}"


class __EmptyAlignmentSegment(AlignmentSegment):
    def __init__(self):
        super().__init__([], 0., Peak.null, [])

    @property
    def startPosition(self):
        return AlignedPair.null

    @property
    def endPosition(self):
        return AlignedPair.null

    def checkForConflicts(self, other: AlignmentSegment):
        return _SegmentPairWithNoConflict(self, other)


AlignmentSegment.empty = __EmptyAlignmentSegment()


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

    def findOptimaPlace(self, conf1, conf2) -> Tuple[AlignmentSegment, AlignmentSegment]:
        if len(conf1[0]) == len(conf2[0]):
            points_1_np = np.cumsum([0] + conf1[1])
            # @todo
            # Check if they should be reversed
            points_2_np = np.cumsum([0] + conf2[1][::-1])[::-1]
            scores_sum = np.add(points_1_np, points_2_np)
            max_index = np.argmax(scores_sum)
            if max_index == 0:
                return self.segment1 - self.conflictingSubsegment1, self.segment2
            elif max_index == len(conf1[0]):
                return self.segment1, self.segment2 - self.conflictingSubsegment2
            else:
                new_seg1 = self.segment1 - \
                    AlignmentSegment(self.conflictingSubsegment1.positions[conf1[2][max_index]:],
                                     sum(p.score for p in
                                         self.conflictingSubsegment1.positions[conf1[2][max_index]:]),
                                     self.segment1.peak)
                new_seg2 = self.segment2 - \
                    AlignmentSegment(self.conflictingSubsegment2.positions[:conf2[2][max_index]],
                                     sum(p.score for p in
                                         self.conflictingSubsegment2.positions[:conf2[2][max_index]]),
                                     self.segment2.peak)
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
            return self.findOptimaPlace(self.conflictingSubsegment1.get_reference_labels(),
                                        self.conflictingSubsegment2.get_reference_labels())
        else:
            return self.findOptimaPlace(self.conflictingSubsegment1.get_query_labels(),
                                        self.conflictingSubsegment2.get_query_labels())
