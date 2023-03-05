from __future__ import annotations

import sys

from p_tqdm import p_imap

from src.alignment.aligner import Aligner, AlignerEngine
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResults
from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.args import Args
from src.correlation.optical_map import OpticalMap
from src.correlation.sequence_generator import SequenceGenerator
from src.diagnostic.diagnostics import DiagnosticsWriter, PrimaryCorrelationDiagnosticsHandler, \
    SecondaryCorrelationDiagnosticsHandler, AlignmentPlotter
from src.messaging.dispatcher import Dispatcher
from src.messaging.messages import InitialAlignmentMessage, CorrelationResultMessage, AlignmentResultRowMessage
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairWithDistanceParser
from src.parsers.xmap_reader import XmapReader


def main():
    args = Args.parse()
    Program(args).run()


class Program:
    def __init__(self, args: Args):
        self.args = args
        self.__readMaps()
        self.xmapReader = XmapReader(XmapAlignmentPairWithDistanceParser(self.referenceMaps, self.queryMaps))
        self.primaryGenerator = SequenceGenerator(args.primaryResolution, args.primaryBlur)
        self.secondaryGenerator = SequenceGenerator(args.secondaryResolution, args.secondaryBlur)
        scorer = AlignmentPositionScorer(args.perfectMatchScore, args.scoreMultiplier, args.unmatchedPenalty)
        segmentsFactory = AlignmentSegmentsFactory(args.minScore, args.breakSegmentThreshold)
        alignerEngine = AlignerEngine(args.maxDistance)
        alignmentSegmentConflictResolver = AlignmentSegmentConflictResolver(SegmentChainer())
        self.aligner = Aligner(scorer, segmentsFactory, alignerEngine, alignmentSegmentConflictResolver)
        self.dispatcher = Dispatcher([])
        if args.diagnosticsEnabled:
            writer = DiagnosticsWriter(args.outputFile)
            self.dispatcher.addHandler(PrimaryCorrelationDiagnosticsHandler(writer))
            self.dispatcher.addHandler(SecondaryCorrelationDiagnosticsHandler(writer))
            self.dispatcher.addHandler(AlignmentPlotter(writer, self.xmapReader, args.benchmarkAlignmentFile))

    def run(self):
        alignmentResultRows = [a for a in p_imap(lambda x: self.__align(*x),
                                                 list((r, q) for r in self.referenceMaps for q in self.queryMaps),
                                                 num_cpus=self.args.numberOfCpus) if a is not None and a.alignedPairs]

        alignmentResult = AlignmentResults.create(self.args.referenceFile.name, self.args.queryFile.name,
                                                  alignmentResultRows)
        self.xmapReader.writeAlignments(self.args.outputFile, alignmentResult)
        if self.args.outputFile is not sys.stdout:
            self.args.outputFile.close()

    def __align(self, referenceMap: OpticalMap, queryMap: OpticalMap):
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator)
        primaryCorrelationReverse = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                                 reverseStrand=True)
        bestPrimaryCorrelation = sorted([primaryCorrelation, primaryCorrelationReverse], key=lambda c: c.getScore())[-1]
        self.dispatcher.dispatch(InitialAlignmentMessage(bestPrimaryCorrelation))

        if not any(bestPrimaryCorrelation.peaks):
            return None

        secondaryCorrelation = bestPrimaryCorrelation.refine(self.secondaryGenerator, self.args.adjustment,
                                                             self.args.peakHeightThreshold)
        self.dispatcher.dispatch(CorrelationResultMessage(bestPrimaryCorrelation, secondaryCorrelation))

        alignmentResultRow = self.aligner.align(referenceMap, queryMap, secondaryCorrelation.peaks,
                                                secondaryCorrelation.reverseStrand)
        self.dispatcher.dispatch(
            AlignmentResultRowMessage(referenceMap, queryMap, alignmentResultRow, bestPrimaryCorrelation))
        return alignmentResultRow

    def __readMaps(self):
        cmapReader = CmapReader()
        with self.args.referenceFile:
            self.referenceMaps = cmapReader.readReferences(self.args.referenceFile, self.args.referenceIds)
        with self.args.queryFile:
            self.queryMaps = cmapReader.readQueries(self.args.queryFile, self.args.queryIds)


if __name__ == '__main__':
    main()
