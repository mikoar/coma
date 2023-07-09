from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List

from p_tqdm import p_imap

from src.alignment.aligner import Aligner
from src.alignment.alignment_results import AlignmentResultRow
from src.args import Args
from src.correlation.optical_map import OpticalMap
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


class SimpleApplicationService(ApplicationService):
    def execute(self, referenceMaps: List[OpticalMap], queryMaps: List[OpticalMap]) -> List[AlignmentResultRow]:
        return [a for a in p_imap(
            lambda x: self.__align(*x),
            list((r, q) for r in referenceMaps for q in queryMaps),
            num_cpus=self.args.numberOfCpus)
                if a is not None and a.alignedPairs]

    def __align(self, referenceMap: OpticalMap, queryMap: OpticalMap) -> AlignmentResultRow | None:
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator)
        primaryCorrelationReverse = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                                 reverseStrand=True)
        bestPrimaryCorrelation = sorted([primaryCorrelation, primaryCorrelationReverse], key=lambda c: c.getScore())[-1]
        self.dispatcher.dispatch(InitialAlignmentMessage(bestPrimaryCorrelation))

        if not any(bestPrimaryCorrelation.peaks):
            return None

        secondaryCorrelation = bestPrimaryCorrelation.refine(self.secondaryGenerator, self.args.secondaryMargin,
                                                             self.args.peakHeightThreshold)
        self.dispatcher.dispatch(CorrelationResultMessage(bestPrimaryCorrelation, secondaryCorrelation))

        alignmentResultRow = self.aligner.align(referenceMap, queryMap, secondaryCorrelation.peaks,
                                                secondaryCorrelation.reverseStrand)
        self.dispatcher.dispatch(
            AlignmentResultRowMessage(referenceMap, queryMap, alignmentResultRow, bestPrimaryCorrelation))
        return alignmentResultRow
