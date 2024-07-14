from __future__ import annotations

from itertools import chain
from typing import List, Iterator

from p_tqdm import p_imap

from multi_pass_workflow_coordinator import MultiPassWorkflowCoordinator
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


class WorkflowCoordinator:
    def __init__(self, args: Args, primaryGenerator: SequenceGenerator, secondaryGenerator: SequenceGenerator,
                 aligner: Aligner, dispatcher: Dispatcher, peaksSelector: PeaksSelector):
        self.args = args
        self.primaryGenerator = primaryGenerator
        self.secondaryGenerator = secondaryGenerator
        self.aligner = aligner
        self.dispatcher = dispatcher
        self.peaksSelector = peaksSelector

    @staticmethod
    def create(args: Args, dispatcher: Dispatcher):
        primaryGenerator = SequenceGenerator(args.primaryResolution, args.primaryBlur)
        secondaryGenerator = SequenceGenerator(args.secondaryResolution, args.secondaryBlur)
        scorer = AlignmentPositionScorer(args.perfectMatchScore, args.distancePenaltyMultiplier, args.unmatchedPenalty)
        segmentsFactory = AlignmentSegmentsFactory(args.minScore, args.breakSegmentThreshold)
        alignerEngine = AlignerEngine(args.maxPairDistance)
        alignmentSegmentConflictResolver = AlignmentSegmentConflictResolver(SegmentChainer())
        aligner = Aligner(scorer, segmentsFactory, alignerEngine, alignmentSegmentConflictResolver)
        if args.outputMode == "single":
            return WorkflowCoordinator(
                args, primaryGenerator, secondaryGenerator, aligner, dispatcher, PeaksSelector(args.peaksCount))
        else:
            return MultiPassWorkflowCoordinator(
                args, primaryGenerator, secondaryGenerator, aligner, dispatcher, PeaksSelector(args.peaksCount))

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
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                          self.args.minPeakDistance, self.args.peaksCount)
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

    def __getAlignmentRow(self, ic: InitialAlignment, sc: CorrelationResult, index: int):
        alignmentResultRow = self.aligner.align(sc.reference, sc.query, sc.peaks, sc.reverseStrand)
        message = AlignmentResultRowMessage(sc.reference, sc.query, alignmentResultRow, ic, index)
        self.dispatcher.dispatch(message)
        return alignmentResultRow, message

    @staticmethod
    def __getBestAlignment(alignmentResultRows: List[AlignmentResultRow]):
        return next(iter(sorted(alignmentResultRows, key=lambda a: a.confidence, reverse=True)), None)
