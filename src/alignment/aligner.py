from __future__ import annotations

import itertools
from itertools import dropwhile, takewhile, chain
from typing import List, NamedTuple

from src.alignment.chained_alignment_position import MockPosition, ChainedAlignmentPosition
from src.alignment.alignment_position import AlignedPair, NotAlignedQueryPosition, \
    NotAlignedReferencePosition, NotAlignedPosition, AlignmentPosition, ScoredAlignedPair, ScoredNotAlignedPosition
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.chain_scorer import ChainScorer
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver, AlignmentSegmentsWithResolvedConflicts
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import OpticalMap, PositionWithSiteId
from src.correlation.peak import Peak


class _ReferenceIndexWithDistance(NamedTuple):
    index: int
    distance: int
    queryShift: int

    @staticmethod
    def withQueryAfterReference(index: int, distance: int):
        return _ReferenceIndexWithDistance(index, distance, distance)

    @staticmethod
    def withQueryBeforeReference(index: int, distance: int):
        return _ReferenceIndexWithDistance(index, distance, -distance)


class AlignerEngine:
    def __init__(self, maxDistance: int):
        self.maxDistance = maxDistance
        self.iteration = 1

    def align(self, reference: OpticalMap, query: OpticalMap, referenceStartPosition: int, referenceEndPosition: int,
              isReverse: bool) -> List[AlignmentPosition]:
        referencePositions = self.__getReferencePositionsWithinRange(reference, referenceStartPosition,
                                                                     referenceEndPosition)
        queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        # alignedPairs = self.__getAlignedPairs(referencePositions, queryPositions, referenceStartPosition)
        # deduplicatedAlignedPairs = list(AlignedPair.deduplicate(alignedPairs))
        # notAlignedPositions = self.__getNotAlignedPositions(queryPositions, referencePositions,
        #                                                     deduplicatedAlignedPairs, referenceStartPosition)
        # return sorted(chain(deduplicatedAlignedPairs, notAlignedPositions))
        return self.__getAlignment(referencePositions, queryPositions, referenceStartPosition)

    def __getReferencePositionsWithinRange(self, reference: OpticalMap, referenceStartPosition: int,
                                           referenceEndPosition: int):
        return list(takewhile(lambda x: x.position <= referenceEndPosition + self.maxDistance, dropwhile(
            lambda x: x.position < referenceStartPosition - self.maxDistance,
            reference.getPositionsWithSiteIds())))

    def __getAlignedPairs(self, referencePositions: List[PositionWithSiteId],
                          queryPositions: List[PositionWithSiteId], referenceStartPosition: int):
        for referencePosition in referencePositions:
            referencePositionAdjustedToQuery = referencePosition.position - referenceStartPosition
            queryPositionsWithinDistance = takewhile(
                lambda x: x.position <= referencePositionAdjustedToQuery + self.maxDistance, dropwhile(
                    lambda x: x.position < referencePositionAdjustedToQuery - self.maxDistance,
                    queryPositions))
            for queryPosition in queryPositionsWithinDistance:
                yield AlignedPair(referencePosition, queryPosition,
                                  queryPosition.position - referencePositionAdjustedToQuery, self.iteration)
        self.iteration += 1
    def __getAlignment(self, referencePositions: List[PositionWithSiteId],
                          queryPositions: List[PositionWithSiteId], referenceStartPosition: int):
        referNAPositions = [NotAlignedReferencePosition(rPos) for rPos in referencePositions]
        queryNAPositions = [NotAlignedQueryPosition(rPos, referenceStartPosition) for rPos in queryPositions]
        allNAPositions = sorted(referNAPositions+queryNAPositions)
        NAScore = -self.maxDistance/6
        ALScore = 2*self.maxDistance/3
        XXScore = -self.maxDistance
        def pairScore(lastPosition, currPosition):
            if isinstance(lastPosition, NotAlignedQueryPosition) == isinstance(currPosition, NotAlignedQueryPosition):
                return XXScore
            else:
                return ALScore-(currPosition.absolutePosition-lastPosition.absolutePosition)
        cumScore = [0, NAScore]
        lastPosition = allNAPositions[0]
        for currPosition in allNAPositions[1:]:
            cumScore.append(max(cumScore[-1]+NAScore, cumScore[-2]+pairScore(lastPosition, currPosition)))
            lastPosition = currPosition
        revAlignmentPositions = []
        i = len(allNAPositions)
        while i>0:
            if cumScore[i] > cumScore[i-1]+NAScore:
                if isinstance(allNAPositions[i-1], NotAlignedQueryPosition):
                    queryPos = allNAPositions[i-1]
                    referPos = allNAPositions[i-2]
                else:
                    queryPos = allNAPositions[i-2]
                    referPos = allNAPositions[i-1]
                queryShift = queryPos.absolutePosition-referPos.absolutePosition
                revAlignmentPositions.append(AlignedPair(referPos.reference, queryPos.query, queryShift))
                i -= 2
            else:
                revAlignmentPositions.append(allNAPositions[i-1])
                i -= 1
        return reversed(revAlignmentPositions)

    @staticmethod
    def __getNotAlignedPositions(queryPositions: List[PositionWithSiteId],
                                 referencePositions: List[PositionWithSiteId],
                                 alignedPairs: List[AlignedPair],
                                 referenceStartPosition: int):
        alignedReferenceSiteIds = [p.reference.siteId for p in alignedPairs]
        alignedQuerySiteIds = [p.query.siteId for p in alignedPairs]
        notAlignedReferencePositions: List[NotAlignedPosition] = \
            [NotAlignedReferencePosition(r) for r in referencePositions if
             r.siteId not in alignedReferenceSiteIds]
        notAlignedQueryPositions = [NotAlignedQueryPosition(q, referenceStartPosition) for q in queryPositions if
                                    q.siteId not in alignedQuerySiteIds]
        return notAlignedReferencePositions + notAlignedQueryPositions


class Aligner:
    def __init__(self, scorer: AlignmentPositionScorer,
                 segmentsFactory: AlignmentSegmentsFactory,
                 alignmentEngine: AlignerEngine,
                 segmentConflictResolver: AlignmentSegmentConflictResolver) -> None:
        self.scorer = scorer
        self.segmentsFactory = segmentsFactory
        self.alignmentEngine = alignmentEngine
        self.segmentConflictResolver = segmentConflictResolver

    def align(self, reference: OpticalMap, query: OpticalMap, peaks: Peak | List[Peak],
              isReverse: bool = False) -> AlignmentResultRow:
        if isinstance(peaks, Peak):
            peaks = [peaks]
        segments = list(itertools.chain.from_iterable(
            [self.getSegments(isReverse, p, query, reference) for p in peaks]))
        queryPositions = list(query.getPositionsWithSiteIds(isReverse))

        return AlignmentResultRow.create(self.segmentConflictResolver.resolveConflicts(segments, queryPositions),
                                         query.moleculeId,
                                         reference.moleculeId,
                                         query.length,
                                         reference.length,
                                         isReverse)

    def getSegments(self, isReverse: bool, peak: Peak, query: OpticalMap, reference: OpticalMap):
        referenceStartPosition = peak.position
        referenceEndPosition = peak.position + query.length
        alignmentPositions = self.alignmentEngine.align(reference, query, referenceStartPosition,
                                                        referenceEndPosition, isReverse)
        scoredPositions = self.scorer.getScoredPositions(alignmentPositions)
        return self.segmentsFactory.getSegments(scoredPositions, peak)




class PeakAligner:
    def __init__(self, scorer: AlignmentPositionScorer, reference: OpticalMap, query: OpticalMap, isReverse: bool, minPeak: int, maxPeak: int):
        self.scorer = scorer
        self.maxDistance = (scorer.perfectMatchScore-2*scorer.unmatchedPenalty)/scorer.distancePenaltyMultiplier
        minRefPos = minPeak - self.maxDistance
        maxRefPos = maxPeak + self.maxDistance + query.length
        self.referencePositions = self.__getReferencePositionsWithinRange(reference, minRefPos, maxRefPos)
        self.queryPositions = list(query.getPositionsWithSiteIds(isReverse))
        self.isReverse = isReverse
        self.queryLength = query.length

    def __getReferencePositionsWithinRange(self, reference: OpticalMap, minRefPos: int, maxRefPos: int):
        return list(takewhile(lambda x: x.position <= maxRefPos, 
                      dropwhile(lambda x: x.position < minRefPos, reference.getPositionsWithSiteIds())))

    def align(self, referenceStartPosition: int) -> (List[AlignmentPosition], List[float]):
        referNAPositions = [ChainedAlignmentPosition(NotAlignedReferencePosition(rPos), \
                                                     self.scorer.unmatchedPenalty, \
                                                     referenceStartPosition) \
                            for rPos in self.referencePositions]
        queryNAPositions = [ChainedAlignmentPosition(NotAlignedQueryPosition(qPos, referenceStartPosition), \
                                                     self.scorer.unmatchedPenalty, \
                                                     referenceStartPosition) \
                            for qPos in self.queryPositions]
        allNAPositions = sorted(referNAPositions+queryNAPositions)
        NAScore = self.scorer.unmatchedPenalty
        def pairScore(lastPosition, currPosition, scorer=self.scorer):
            if lastPosition.isNAQuery() == currPosition.isNAQuery():
                return 3*scorer.unmatchedPenalty
            else:
                return scorer.perfectMatchScore - scorer.distancePenaltyMultiplier*(currPosition.queryPosition-lastPosition.queryPosition)
        cumScoreAll = [0, NAScore]
        lastPosition = allNAPositions[0]
        for currPosition in allNAPositions[1:]:
            cumScoreAll.append(max(cumScoreAll[-1]+NAScore, cumScoreAll[-2]+pairScore(lastPosition, currPosition)))
            lastPosition = currPosition
        revAlignmentPositions = []
        i = len(allNAPositions)
        while i>0:
            if cumScoreAll[i] > cumScoreAll[i-1]+NAScore:
                if allNAPositions[i-1].isNAQuery():
                    queryPos = allNAPositions[i-1]
                    referPos = allNAPositions[i-2]
                else:
                    queryPos = allNAPositions[i-2]
                    referPos = allNAPositions[i-1]
                posShift = queryPos.queryPosition-referPos.queryPosition
                ap = AlignedPair(referPos.alignmentPosition.reference, queryPos.alignmentPosition.query, posShift)
                score = cumScoreAll[i] - cumScoreAll[i-2]
                revAlignmentPositions.append(ChainedAlignmentPosition(ap, score, referenceStartPosition))
                i -= 2
            else:
                revAlignmentPositions.append(allNAPositions[i-1])
                i -= 1
        for i, ap in enumerate(revAlignmentPositions):
            if ap.queryPosition<=self.queryLength-self.maxDistance or not ap.isNARefer():
                break
        assert all(ap.isNARefer() for ap in revAlignmentPositions[:i])
        alignmentPositions = list(reversed(revAlignmentPositions[i:]))
        for i, ap in enumerate(alignmentPositions):
            if ap.queryPosition>self.maxDistance or not ap.isNARefer():
                break
        assert all(ap.isNARefer() for ap in alignmentPositions[:i])
        startPosition = ChainedAlignmentPosition(MockPosition(0, True), 
                                                 self.scorer.endReachingScore, referenceStartPosition)
        endPosition = ChainedAlignmentPosition(MockPosition(self.queryLength+1, False), 
                                                 self.scorer.endReachingScore, referenceStartPosition)
        # startPosition = ChainedAlignmentPosition(MockPosition(alignmentPositions[i].queryPosition-1, True), 
        #                                          self.scorer.endReachingScore, referenceStartPosition)
        # endPosition = ChainedAlignmentPosition(MockPosition(alignmentPositions[-1].queryPosition+1, False), 
        #                                          self.scorer.endReachingScore, referenceStartPosition)
        return [startPosition] + alignmentPositions[i:] + [endPosition]


class ChainBuilder:
    def __init__(self, alPosScorer: AlignmentPositionScorer, chainScorer: ChainScorer, blurredEnd: int) -> None:
        self.alPosScorer = alPosScorer
        self.chainScorer = chainScorer
        self.blurredEnd = blurredEnd+1

    def setFirstPass(self, isFirstPass):
        self.chainScorer.setFirstPass(isFirstPass)

    def align(self, reference: OpticalMap, query: OpticalMap, peaks: List[Peak],
              isReverse: bool = False) -> AlignmentResultRow:
        if not peaks:
            # return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts([]),
            #                                  query.moleculeId,
            #                                  reference.moleculeId,
            #                                  query.length,
            #                                  reference.length,
            #                                  isReverse)
            return AlignmentResultRow([], 
                                      query.moleculeId,
                                      reference.moleculeId,
                                      query.length,
                                      reference.length,
                                      reverseStrand = isReverse)
        peakShifts = [peak.position for peak in peaks]
        peakAligner = PeakAligner(self.alPosScorer, reference, query, isReverse, min(peakShifts), max(peakShifts))
        PeakAlignmentPositions = {}
        AllAlignmentPositions = []
        for shift in peakShifts:
            alignmentPositions = peakAligner.align(shift)
            PeakAlignmentPositions[shift] = alignmentPositions
            AllAlignmentPositions.extend(alignmentPositions)
        AllAlignmentPositions.sort()

        peakAPIndexes = {shift: [] for shift in peakShifts}
        bestPrevIndexes = []
        AllCumScores = []
        for i, currAP in enumerate(AllAlignmentPositions):
            currShift = currAP.queryShift
            if currAP.isStart():
                assert not peakAPIndexes[currShift]
                peakAPIndexes[currShift] = [i]
                bestPrevIndexes.append(-1)
                AllCumScores.append(currAP.score)
                continue
            assert peakAPIndexes[currShift]
            if currAP.isEnd():
                prev = peakAPIndexes[currShift][-1]
                peakAPIndexes[currShift].append(i)
                bestPrevIndexes.append(prev)
                AllCumScores.append(AllCumScores[prev] + currAP.score)
                continue
            currScores = [(0, -1)]
            currShiftLastAP = AllAlignmentPositions[peakAPIndexes[currShift][-1]]
            for prevShift in peakShifts:
                if prevShift==currShift:
                    prev = peakAPIndexes[prevShift][-1]
                    currScores.append((AllCumScores[prev] + currAP.score, prev))
                else:
                    prevCandidates = []
                    for prevCand in peakAPIndexes[prevShift][::-1]:
                        prevCandAP = AllAlignmentPositions[prevCand]
                        if prevCandAP.isStart() or prevCandAP.lessOnBothMaps(currShiftLastAP):
                            break
                        if prevCandAP.lessOnBothMaps(currAP):
                            prevCandidates.append((AllCumScores[prevCand], prevCand))
                    if prevCandidates:
                        prevScore, prev = max(prevCandidates)
                        indelLength = abs(prevShift-currShift)
                        indelScore = self.chainScorer.scoreIndel(indelLength)
                        currScores.append((prevScore + indelScore + currAP.score, prev))
            peakAPIndexes[currShift].append(i)
            bestScore, bestPrev = max(currScores)
            AllCumScores.append(bestScore)
            bestPrevIndexes.append(bestPrev)
            
        maxScore = 0
        maxIndex = None
        for i, cumScore in enumerate(AllCumScores):
            if cumScore > maxScore:
                maxScore = cumScore
                maxIndex = i
        assert maxScore>0 and maxIndex is not None
        if maxScore < self.chainScorer.minAlignmentScore():
            # return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts([]),
            #                                  query.moleculeId,
            #                                  reference.moleculeId,
            #                                  query.length,
            #                                  reference.length,
            #                                  isReverse)
            return AlignmentResultRow([], 
                                      query.moleculeId,
                                      reference.moleculeId,
                                      query.length,
                                      reference.length,
                                      reverseStrand = isReverse)
           
        currAP = AllAlignmentPositions[maxIndex]
        if currAP.isStart():
            # print("Warning: no positive alignment score for", query.moleculeId, isReverse, query.length, peakShifts)
            # return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts([]),
            #                                  query.moleculeId,
            #                                  reference.moleculeId,
            #                                  query.length,
            #                                  reference.length,
            #                                  isReverse)
            return AlignmentResultRow([], 
                                      query.moleculeId,
                                      reference.moleculeId,
                                      query.length,
                                      reference.length,
                                      reverseStrand = isReverse)
        lastAP = currAP
        revChain = []
        # revIndels = []
        rests = []
        revSegment = []
        if not lastAP.isEnd():
            revSegment.append(lastAP)
            alignmentEnd = max(lastAP.queryPosition+1, AllAlignmentPositions[maxIndex+1].queryPosition - self.blurredEnd)
            rests.append(query.getSubMap(isReverse, start = alignmentEnd))
            # rests.append(list(filter(lambda pos: pos.position>lastAP.queryPosition, query.getPositionsWithSiteIds(isReverse))))
        # revSegment = [] if currAP.isEnd() else [currAP]
        segmentShift = currAP.queryShift
        currIndex = bestPrevIndexes[maxIndex]
        assert not (currAP.isEnd() and AllAlignmentPositions[currIndex].queryShift != segmentShift)
        while currIndex>-1:
            currAP = AllAlignmentPositions[currIndex]
            assert not currAP.isEnd()
            if currAP.queryShift == segmentShift:
                if not currAP.isStart():
                    revSegment.append(currAP)
            else:
                assert revSegment
                revChain.append((segmentShift, reversed(revSegment)))
                # revIndels.append((currAP, revSegment[-1]))
                assert not currAP.isStart()
                revSegment = [currAP]
                segmentShift = currAP.queryShift
            maxIndex = currIndex
            currIndex = bestPrevIndexes[currIndex]
        assert revSegment
        revChain.append((segmentShift, reversed(revSegment)))
        firstAP = currAP
        if not firstAP.isStart():
            alignmentStart = min(firstAP.queryPosition-1, AllAlignmentPositions[maxIndex-1].queryPosition + self.blurredEnd)
            rests.append(query.getSubMap(isReverse, end = alignmentStart))
            # rests.append(list(filter(lambda pos: pos.position<firstAP.queryPosition, query.getPositionsWithSiteIds(isReverse))))
        
        def smapData(ap):
            refPos = ap.reference.position
            queryPos = query.getAbsolutePosition(ap.query.position, isReverse)
            refIdx = ap.reference.siteId
            queryIdx = ap.query.siteId
            return refPos, queryPos, refIdx, queryIdx
            
        
        shiftToPeak = {peak.position: peak for peak in peaks}
        finalSegments = []
        indels = []
        prevAP = None
        for shift, segment in reversed(revChain):
            positions = [ScoredAlignedPair(ap.alignmentPosition, ap.score) if ap.isAlPair() else
                         ScoredNotAlignedPosition(ap.alignmentPosition, ap.score) for ap in segment]
            APs = [ap for ap in positions if isinstance(ap, AlignedPair)]
            if APs:
                if prevAP is not None:
                    indels.append((smapData(prevAP), smapData(APs[0]), shift-prevShift))
                prevAP = APs[-1]
                prevShift = shift
            peak = shiftToPeak[segmentShift]
            allPeakPositions = []
            finalSegments.append(AlignmentSegment.create(positions, peak, allPeakPositions))


        # return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts(finalSegments),
        #                                  query.moleculeId,
        #                                  reference.moleculeId,
        #                                  query.length,
        #                                  reference.length,
        #                                  isReverse)

        queryStartPosition = query.getAbsolutePosition(firstAP.queryPosition, isReverse)
        queryEndPosition =  query.getAbsolutePosition(lastAP.queryPosition, isReverse)
        if isReverse:
            queryStartPosition, queryEndPosition = queryEndPosition, queryStartPosition
        referenceStartPosition = firstAP.referencePosition()
        referenceEndPosition = lastAP.referencePosition()
        return AlignmentResultRow(finalSegments, 
                                  query.moleculeId,
                                  reference.moleculeId,
                                  query.length,
                                  reference.length,
                                  queryStartPosition,
                                  queryEndPosition, 
                                  referenceStartPosition, 
                                  referenceEndPosition, 
                                  isReverse,
                                  maxScore, 
                                  indels = indels, 
                                  rests = rests)

