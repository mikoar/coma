from __future__ import annotations

from abc import ABC, abstractmethod
from itertools import chain
from typing import List, Iterator

from p_tqdm import p_imap

from src.alignment.aligner import Aligner, AlignerEngine
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.args import Args
from src.correlation.optical_map import OpticalMap, InitialAlignment, CorrelationResult
from src.correlation.peaks_selector import PeaksSelector, SelectedPeak
from src.correlation.sequence_generator import SequenceGenerator
from src.extensions.dispatcher import Dispatcher
from src.extensions.messages import CorrelationResultMessage, InitialAlignmentMessage, AlignmentResultRowMessage, \
    MultipleAlignmentResultRowsMessage


class ApplicationServiceFactory:
    def create(self, args: Args, dispatcher: Dispatcher) -> ApplicationService:
        primaryGenerator = SequenceGenerator(args.primaryResolution, args.primaryBlur)
        secondaryGenerator = SequenceGenerator(args.secondaryResolution, args.secondaryBlur)
        scorer = AlignmentPositionScorer(args.perfectMatchScore, args.distancePenaltyMultiplier, args.unmatchedPenalty)
        segmentsFactory = AlignmentSegmentsFactory(args.minScore, args.breakSegmentThreshold)
        alignerEngine = AlignerEngine(args.maxDistance)
        alignmentSegmentConflictResolver = AlignmentSegmentConflictResolver(SegmentChainer())
        aligner = Aligner(scorer, segmentsFactory, alignerEngine, alignmentSegmentConflictResolver)
        if args.onePeakPerReference:
            return OnePeakPerReferenceApplicationService(
                args, primaryGenerator, secondaryGenerator, aligner, dispatcher)
        else:
            return MultiPeakApplicationService(
                args, primaryGenerator, secondaryGenerator, aligner, dispatcher, PeaksSelector(args.peaksCount))


class ApplicationService(ABC):
    def __init__(self,
                 args: Args,
                 primaryGenerator: SequenceGenerator,
                 secondaryGenerator: SequenceGenerator,
                 aligner: Aligner,
                 dispatcher: Dispatcher):
        self.args = args
        self.primaryGenerator = primaryGenerator
        self.secondaryGenerator = secondaryGenerator
        self.aligner = aligner
        self.dispatcher = dispatcher

    @abstractmethod
    def execute(self, referenceMaps: List[OpticalMap], queryMaps: List[OpticalMap]) -> List[AlignmentResultRow]:
        pass


class OnePeakPerReferenceApplicationService(ApplicationService):
    def execute(self, referenceMaps: List[OpticalMap], queryMaps: List[OpticalMap]) -> List[AlignmentResultRow]:
        return [a for a in p_imap(
            lambda x: self.__align(*x),
            list((r, q) for r in referenceMaps for q in queryMaps),
            num_cpus=self.args.numberOfCpus,
            disable=self.args.disableProgressBar)
                if a is not None and a.alignedPairs]

    def __align(self, referenceMap: OpticalMap, queryMap: OpticalMap) -> AlignmentResultRow | None:
        primaryCorrelation = self.__getPrimaryCorrelation(referenceMap, queryMap)
        if not any(primaryCorrelation.peaks):
            return None

        secondaryCorrelation = primaryCorrelation.refine(primaryCorrelation.maxPeak.position, self.secondaryGenerator,
                                                         self.args.secondaryMargin, self.args.peakHeightThreshold)
        self.dispatcher.dispatch(CorrelationResultMessage(primaryCorrelation, secondaryCorrelation))

        alignmentResultRow = self.aligner.align(referenceMap, queryMap, secondaryCorrelation.peaks,
                                                secondaryCorrelation.reverseStrand)
        self.dispatcher.dispatch(
            AlignmentResultRowMessage(referenceMap, queryMap, alignmentResultRow, primaryCorrelation))
        return alignmentResultRow

    def __getPrimaryCorrelation(self, referenceMap: OpticalMap, queryMap: OpticalMap):
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                          self.args.minPeakDistance)
        self.dispatcher.dispatch(InitialAlignmentMessage(primaryCorrelation))
        primaryCorrelationReverse = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                                 self.args.minPeakDistance, reverseStrand=True)
        self.dispatcher.dispatch(InitialAlignmentMessage(primaryCorrelationReverse))
        bestPrimaryCorrelation = sorted([primaryCorrelation, primaryCorrelationReverse], key=lambda c: c.getScore())[-1]
        return bestPrimaryCorrelation


class MultiPeakApplicationService(ApplicationService):
    def __init__(self, args: Args, primaryGenerator: SequenceGenerator, secondaryGenerator: SequenceGenerator,
                 aligner: Aligner, dispatcher: Dispatcher, peaksSelector: PeaksSelector):
        super().__init__(args, primaryGenerator, secondaryGenerator, aligner, dispatcher)
        self.peaksSelector = peaksSelector

    def execute(self, referenceMaps: List[OpticalMap], queryMaps: List[OpticalMap]) -> List[AlignmentResultRow]:
        return [a for a in p_imap(
            lambda x: self.__align(*x),
            list((referenceMaps, q) for q in queryMaps),
            num_cpus=self.args.numberOfCpus,
            disable=self.args.disableProgressBar)
                if a is not None and a.alignedPairs]

    def __align(self, referenceMaps: List[OpticalMap], queryMap: OpticalMap) -> AlignmentResultRow | None:
        primaryCorrelations = chain.from_iterable(self.__getPrimaryCorrelations(r, queryMap) for r in referenceMaps)

        bestPrimaryCorrelationPeaks = self.peaksSelector.selectPeaks(primaryCorrelations)

        secondaryCorrelations = [self.__getSecondaryCorrelation(p, i)
                                 for i, p in enumerate(bestPrimaryCorrelationPeaks)]

        alignmentResultRows, messages = zip(*[self.__getAlignmentRow(pc, sc, i) for i, (pc, sc) in
                                              enumerate(secondaryCorrelations)])
        self.dispatcher.dispatch(MultipleAlignmentResultRowsMessage(messages))
        return self.__getBestAlignment(alignmentResultRows)

    def __getPrimaryCorrelations(self, referenceMap: OpticalMap, queryMap: OpticalMap) -> Iterator[InitialAlignment]:
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator, self.args.minPeakDistance, self.args.peaksCount)
        self.dispatcher.dispatch(InitialAlignmentMessage(primaryCorrelation))

        primaryCorrelationReverse = queryMap.getInitialAlignment(
            referenceMap, self.primaryGenerator, self.args.minPeakDistance, self.args.peaksCount, reverseStrand=True)

        self.dispatcher.dispatch(InitialAlignmentMessage(primaryCorrelationReverse))
        if any(primaryCorrelation.peaks):
            yield primaryCorrelation
        if any(primaryCorrelationReverse.peaks):
            yield primaryCorrelationReverse

    def __getSecondaryCorrelation(self, selectedPeak: SelectedPeak, index: int):
        secondaryCorrelation = selectedPeak.primaryCorrelation.refine(selectedPeak.peak.position,
                                                                      self.secondaryGenerator,
                                                                      self.args.secondaryMargin,
                                                                      self.args.peakHeightThreshold)

        self.dispatcher.dispatch(CorrelationResultMessage(selectedPeak.primaryCorrelation, secondaryCorrelation, index))
        return selectedPeak.primaryCorrelation, secondaryCorrelation

    def __getAlignmentRow(self, sc: CorrelationResult, ic: InitialAlignment, index: int):
        alignmentResultRow = self.aligner.align(sc.reference, sc.query, sc.peaks, sc.reverseStrand)
        message = AlignmentResultRowMessage(sc.reference, sc.query, alignmentResultRow, ic, index)
        return alignmentResultRow, message

    @staticmethod
    def __getBestAlignment(alignmentResultRows: List[AlignmentResultRow]):
        return next(iter(sorted(alignmentResultRows, key=lambda a: a.confidence, reverse=True)), None)
