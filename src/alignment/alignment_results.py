from __future__ import annotations

import itertools
import argparse
import os
from dataclasses import dataclass
from enum import Enum
from typing import List

from src.alignment.alignment_position import AlignedPair, NotAlignedPosition
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentsWithResolvedConflicts
from src.alignment.segments import AlignmentSegment
from src.diagnostic.benchmark_alignment import BenchmarkAlignment
from src.correlation.optical_map import OpticalMap


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"


@dataclass
class AlignmentResults:
    referenceFilePath: str
    queryFilePath: str
    rows: List[AlignmentResultRow]

    @staticmethod
    def create(referenceFilePath: str,
               queryFilePath: str,
               rows: List[AlignmentResultRow],
               rows_rest: List[AlignmentResultRow],
               mode: str,
               out_file: argparse.FileType,
               maxDifference: int):
        if mode == "best":
            rows = rows + rows_rest
        rowsSortedByQueryIdThenByConfidence = \
            sorted(sorted(rows, key=lambda r: r.confidence, reverse=True), key=lambda r: r.queryId)
        rowsWithoutSubsequentAlignmentsForSingleQuery = \
            [next(group) for _, group in itertools.groupby(rowsSortedByQueryIdThenByConfidence, lambda r: r.queryId)]
        new_file = open("{0}_{2}{1}".format(*os.path.splitext(out_file.name) + (1,)),
                        mode='w', encoding=out_file.encoding)
        rowsSortedByQueryIdThenByConfidenceRest = \
            sorted(sorted(rows_rest, key=lambda r: r.confidence, reverse=True), key=lambda r: r.queryId)
        rowsWithoutSubsequentAlignmentsForSingleQueryRest = \
            [next(group) for _, group in itertools.groupby(rowsSortedByQueryIdThenByConfidenceRest, lambda r: r.queryId)]
        if mode == 'separate':
            return [(out_file, AlignmentResults(referenceFilePath, queryFilePath, rowsWithoutSubsequentAlignmentsForSingleQuery)),
                    (new_file, AlignmentResults(referenceFilePath, queryFilePath, rowsWithoutSubsequentAlignmentsForSingleQueryRest))]
        
        joinedRows, separateRows = AlignmentResults.resolve(rowsWithoutSubsequentAlignmentsForSingleQuery,
                                                            rowsWithoutSubsequentAlignmentsForSingleQueryRest,
                                                            maxDifference)
        if mode == 'best':
            rowsWithoutSubsequentAlignmentsForSingleQueryBest = sorted(joinedRows + rowsWithoutSubsequentAlignmentsForSingleQuery,
                                                                       key=lambda r: r.queryId)
            return [(out_file, AlignmentResults(referenceFilePath, queryFilePath, rowsWithoutSubsequentAlignmentsForSingleQueryBest))]
        if mode == 'joined':
            return [(out_file, AlignmentResults(referenceFilePath, queryFilePath, joinedRows)),
                    (new_file, AlignmentResults(referenceFilePath, queryFilePath,  separateRows))]
        if mode == 'all':
            third_file = open("{0}_{2}{1}".format(*os.path.splitext(out_file.name) + (2,)),
                              mode='w', encoding=out_file.encoding)
            return [(out_file, AlignmentResults(referenceFilePath, queryFilePath, joinedRows)),
                    (new_file, AlignmentResults(referenceFilePath, queryFilePath,  rowsWithoutSubsequentAlignmentsForSingleQuery)),
                    (third_file, AlignmentResults(referenceFilePath, queryFilePath,  rowsWithoutSubsequentAlignmentsForSingleQueryRest))]
        

    @staticmethod
    def resolve(rows: List[AlignmentResultRow],
                rows_rest: List[AlignmentResultRow],
                maxDifference: int):
        separate = []
        joined = []
        for _, queries in itertools.groupby(sorted(rows + rows_rest, key=lambda r: r.referenceId), lambda r: r.referenceId):
            for _, group in itertools.groupby(sorted(list(queries), key=lambda r: r.queryId), lambda r: r.queryId):
                group = list(group)
                if len(group) == 1:
                    separate.append(group[0])
                else:
                    if group[0].check_overlap(group[1], maxDifference):
                        resolved = group[0].resolve(group[1])
                        if resolved:
                            joined.append(resolved)
                        else:
                            separate.extend(group)
                    else:
                        separate.extend(group)
        return joined, separate

class AlignmentResultRow(BenchmarkAlignment):
    @staticmethod
    def create(segmentsWithoutConflicts: AlignmentSegmentsWithResolvedConflicts,
               queryId: int,
               referenceId: int,
               queryLength: int,
               referenceLength: int,
               reverseStrand: bool):
        segments = segmentsWithoutConflicts.segments
        alignedPairs = sorted(p for s in segments for p in s.positions if isinstance(p, AlignedPair))
        firstPair = alignedPairs[0] if alignedPairs else AlignedPair.null
        lastPair = alignedPairs[-1] if alignedPairs else AlignedPair.null
        queryStartPosition = (firstPair if not reverseStrand else lastPair).query.position
        queryEndPosition = (lastPair if not reverseStrand else firstPair).query.position
        referenceStartPosition = firstPair.reference.position
        referenceEndPosition = lastPair.reference.position
        confidence = sum(s.segmentScore for s in segments)
        return AlignmentResultRow(segments, queryId, referenceId, queryLength, referenceLength, queryStartPosition,
                                  queryEndPosition, referenceStartPosition, referenceEndPosition, reverseStrand,
                                  confidence)

    def __init__(self,
                 segments: List[AlignmentSegment],
                 queryId: int = 1,
                 referenceId: int = 1,
                 queryLength: int = 1,
                 referenceLength: int = 1,
                 queryStartPosition: int = 0,
                 queryEndPosition: int = 0,
                 referenceStartPosition: int = 0,
                 referenceEndPosition: int = 0,
                 reverseStrand: bool = False,
                 confidence: float = 0.,
                 alignedRest: bool = False):

        self.queryId = queryId
        self.referenceId = referenceId
        self.queryStartPosition = queryStartPosition
        self.queryEndPosition = queryEndPosition
        self.referenceStartPosition = referenceStartPosition
        self.referenceEndPosition = referenceEndPosition
        self.reverseStrand = reverseStrand
        self.confidence = confidence
        self.queryLength = queryLength
        self.referenceLength = referenceLength
        self.segments = segments
        self.alignedRest = alignedRest

    @property
    def positions(self):
        return [position for segment in self.segments for position in segment.positions]

    @property
    def alignedPairs(self) -> List[AlignedPair]:
        return [p for p in self.positions if isinstance(p, AlignedPair)]

    @property
    def notAlignedPositions(self) -> List[NotAlignedPosition]:
        return [p for p in self.positions if isinstance(p, NotAlignedPosition)]

    @property
    def cigarString(self):
        if not self.alignedPairs:
            return ""
        hitEnums = list(self.__getHitEnums())
        return "".join(self.__aggregateHitEnums(hitEnums))

    def __getHitEnums(self):
        pairs = list(self.__removeDuplicateQueryPositionsPreservingLastOne(self.alignedPairs))
        pairsIterator = iter(pairs)
        currentPair: AlignedPair = next(pairsIterator)
        previousQuery = currentPair.query.siteId
        for referenceIndex in range(pairs[0].reference.siteId,
                                    pairs[-1].reference.siteId + 1):
            queryIncrement = abs(currentPair.query.siteId - previousQuery)
            if queryIncrement > 1:
                for _ in range(1, queryIncrement):
                    yield HitEnum.INSERTION
                previousQuery = currentPair.query.siteId
            if currentPair.reference.siteId == referenceIndex:
                previousQuery = currentPair.query.siteId
                currentPair = next(pairsIterator, None)
                yield HitEnum.MATCH
            elif currentPair.reference.siteId > referenceIndex:
                yield HitEnum.DELETION

    @staticmethod
    def __removeDuplicateQueryPositionsPreservingLastOne(pairs: List[AlignedPair]):
        for _, ambiguousPairs in itertools.groupby(pairs, lambda pair: pair.query.siteId):
            *_, lastPair = ambiguousPairs
            yield lastPair

    @staticmethod
    def __aggregateHitEnums(hits: List[HitEnum]):
        hit = None
        count = 1
        previousHit: HitEnum = hits[0]
        for hit in hits[1:]:
            if hit == previousHit:
                count += 1
            else:
                yield AlignmentResultRow.__hitToString(count, previousHit)
                previousHit = hit
                count = 1
        if hit:
            yield AlignmentResultRow.__hitToString(count, hit)

    @staticmethod
    def __hitToString(count, hit):
        x = f"{count}{hit.value}"
        return x

    def getUnalignedFragments(self, queries:List[OpticalMap]) -> List[OpticalMap]:
        """Function used to return unaligned fragments of the query
        if those parts constitute more than 0.2 of the whole query

        :param queries: whole query which is being currently aligned
        :type queries: List[OpticalMap]
        :return: Unaligned parts of query in question
        :rtype: List[OpticalMap]
        """
        if abs(self.queryStartPosition - self.queryEndPosition) > 0.8  * self.queryLength:
            return []
        else:
            query = next((opticMap for opticMap in queries if opticMap.moleculeId == self.queryId), None)
            if self.queryStartPosition == 0.0 or self.queryEndPosition == 0.0:
                # Aligned positions are at the end/start
                if self.orientation == '+':
                    positions = query.positions[query.positions.index(self.queryEndPosition) - 2 :]
                    shift = len(query.positions)-len(positions)
                    if self.queryEndPosition == 0.0:
                        return []
                else:
                    alignedPairs = sorted(p for s in self.segments for p in s.positions if isinstance(p, AlignedPair))
                    lastPair = alignedPairs[-1] if alignedPairs else AlignedPair.null
                    positions = query.positions[: lastPair.query.siteId + 3]
                    shift = 0
                return [OpticalMap(self.queryId, self.queryLength, positions, shift=shift)]
            else:
                # Case where aligned fragment is in the middle
                if self.orientation == '+':
                    positions1 = query.positions[: query.positions.index(self.queryStartPosition) + 3]
                    positions2 = query.positions[query.positions.index(self.queryEndPosition) - 2 :]

                else:
                    alignedPairs = sorted(p for s in self.segments for p in s.positions if isinstance(p, AlignedPair))
                    firstPair = alignedPairs[0] if alignedPairs else AlignedPair.null
                    lastPair = alignedPairs[-1] if alignedPairs else AlignedPair.null
                    positions1 = query.positions[: lastPair.query.siteId + 3]
                    positions2 = query.positions[firstPair.query.siteId - 2 :]

                if len(positions1) >= 7 and len(positions2) >= 7:
                    return [OpticalMap(self.queryId, self.queryLength, positions1, shift=0),
                            OpticalMap(self.queryId, self.queryLength, positions2, shift=len(query.positions)-len(positions2))]
                elif len(positions1) >= 7:
                    return [OpticalMap(self.queryId, self.queryLength, positions1, shift=0)]
                elif len(positions2) >= 7:
                    return [OpticalMap(self.queryId, self.queryLength, positions2, shift=len(query.positions)-len(positions2))]
                else:
                    return []

    def setAlignedRest(self, alignedRest: bool):
        self.alignedRest = alignedRest
        return self

    def check_overlap(self, alignedRest: AlignmentResultRow, maxDifference:int) -> bool:
        """Function used to identify overlapping alignments of the same query

        :param alignedRest: Other alignment of the same query
        :type alignedRest: AlignmentResultRow
        :param maxDifference: Maximum difference between reference positions of the
        alignments if they are to be joint
        :type maxDifference: int
        :return: Whether those two alignments should be joined
        :rtype: bool
        """
        if self.orientation == alignedRest.orientation and self.referenceId == alignedRest.referenceId:
            diff = abs(max(self.referenceStartPosition, alignedRest.referenceStartPosition) - \
                min(self.referenceEndPosition, alignedRest.referenceEndPosition))
            if diff <= maxDifference:
                return True
        return False

    def resolve(self, alignedRest: AlignmentResultRow) -> AlignmentResultRow:
        """Function used to resolve conflicts between two overlapping alignments of the same query

        :param alignedRest: Other alignment of the same query
        :type alignedRest: AlignmentResultRow
        :return: Joint alignment with resolved conflicts
        :rtype: AlignmentResultRow
        """
        if self.alignedPairs[0].reference.position < alignedRest.alignedPairs[0].reference.position:
            pair = self.segments[0].checkForConflicts(alignedRest.segments[0])
        else:
            pair = alignedRest.segments[0].checkForConflicts(self.segments[0])

        if pair.resolveConflict():
            seg1, seg2 = pair.resolveConflict()
            notEmptySegments = [s for s in [seg1, seg2] if s != AlignmentSegment.empty]
            return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts(notEmptySegments),
                                      self.queryId, self.referenceId, self.queryLength, self.referenceLength,
                                      self.reverseStrand)
        return
