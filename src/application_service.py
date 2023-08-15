from __future__ import annotations

from abc import ABC, abstractmethod
from itertools import chain
from typing import List, Iterator

from p_tqdm import p_imap

from src.alignment.aligner import Aligner
from src.alignment.alignment_results import AlignmentResultRow
from src.args import Args
from src.correlation.optical_map import OpticalMap, InitialAlignment, CorrelationResult
from src.correlation.peaks_selector import PeaksSelector, SelectedPeak
from src.correlation.sequence_generator import SequenceGenerator
from src.messaging.dispatcher import Dispatcher
from src.messaging.messages import CorrelationResultMessage, InitialAlignmentMessage, AlignmentResultRowMessage


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
            num_cpus=self.args.numberOfCpus)
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
            num_cpus=self.args.numberOfCpus)
                if a is not None and a.alignedPairs]

    def __align(self, referenceMaps: List[OpticalMap], queryMap: OpticalMap) -> AlignmentResultRow | None:
        primaryCorrelations = chain.from_iterable(self.__getPrimaryCorrelations(r, queryMap) for r in referenceMaps)

        bestPrimaryCorrelationPeaks = self.peaksSelector.selectPeaks(primaryCorrelations)

        secondaryCorrelations = [self.__getSecondaryCorrelation(p) for p in bestPrimaryCorrelationPeaks]

        alignmentResultRows = [self.__getAlignmentRow(pc, sc) for pc, sc in secondaryCorrelations]
        return self.__getBestAlignment(alignmentResultRows)

    def __getPrimaryCorrelations(self, referenceMap: OpticalMap, queryMap: OpticalMap) -> Iterator[InitialAlignment]:
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                          self.args.minPeakDistance)
        primaryCorrelationReverse = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                                 self.args.minPeakDistance, reverseStrand=True)
        self.dispatcher.dispatch(InitialAlignmentMessage(primaryCorrelation))
        self.dispatcher.dispatch(InitialAlignmentMessage(primaryCorrelationReverse))
        if any(primaryCorrelation.peaks):
            yield primaryCorrelation
        if any(primaryCorrelationReverse.peaks):
            yield primaryCorrelationReverse

    def __getSecondaryCorrelation(self, selectedPeak: SelectedPeak):
        secondaryCorrelation = selectedPeak.primaryCorrelation.refine(selectedPeak.peak.position,
                                                                      self.secondaryGenerator,
                                                                      self.args.secondaryMargin,
                                                                      self.args.peakHeightThreshold)

        self.dispatcher.dispatch(CorrelationResultMessage(selectedPeak.primaryCorrelation, secondaryCorrelation))
        return selectedPeak.primaryCorrelation, secondaryCorrelation

    def __getAlignmentRow(self, sc: CorrelationResult, ic: InitialAlignment):
        alignmentResultRow = self.aligner.align(sc.reference, sc.query, sc.peaks, sc.reverseStrand)
        self.dispatcher.dispatch(AlignmentResultRowMessage(sc.reference, sc.query, alignmentResultRow, ic))
        return alignmentResultRow

    @staticmethod
    def __getBestAlignment(alignmentResultRows):
        return next(sorted(alignmentResultRows, key=lambda a: a.confidence, reverse=True), default=None)
