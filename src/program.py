from __future__ import annotations

import sys
from typing import List

from p_tqdm import p_imap

from src.alignment.aligner import Aligner, AlignerEngine
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResults
from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.args import Args
from src.correlation.optical_map import OpticalMap
from src.correlation.plot import plotRefinedCorrelation, plotCorrelation
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader


def main():
    args = Args.parse()
    Program(args).run()


class Program:
    def __init__(self, args: Args):
        self.args = args
        self.cmapReader = CmapReader()
        self.xmapReader = XmapReader()
        self.primaryGenerator = SequenceGenerator(args.primaryResolution, args.primaryBlur)
        self.secondaryGenerator = SequenceGenerator(args.secondaryResolution, args.secondaryBlur)
        scorer = AlignmentPositionScorer(args.perfectMatchScore, args.scoreMultiplier, args.unmatchedPenalty)
        segmentsFactory = AlignmentSegmentsFactory(args.minScore, args.breakSegmentThreshold)
        alignerEngine = AlignerEngine(args.maxDistance)
        alignmentSegmentConflictResolver = AlignmentSegmentConflictResolver(SegmentChainer())
        self.aligner = Aligner(scorer, segmentsFactory, alignerEngine, alignmentSegmentConflictResolver)

    def run(self):
        referenceMaps: List[OpticalMap]
        queryMaps: List[OpticalMap]
        with self.args.referenceFile:
            referenceMaps = self.cmapReader.readReferences(self.args.referenceFile, self.args.referenceIds)
        with self.args.queryFile:
            queryMaps = self.cmapReader.readQueries(self.args.queryFile, self.args.queryIds)

        alignmentResultRows = [r for r in p_imap(lambda x: self.__align(*x),
                                                 list((r, q) for r in referenceMaps for q in queryMaps),
                                                 num_cpus=self.args.numberOfCpus) if r is not None and r.alignedPairs]
        alignmentResult = AlignmentResults(self.args.referenceFile.name, self.args.queryFile.name, alignmentResultRows)

        self.xmapReader.writeAlignments(self.args.outputFile, alignmentResult)
        if self.args.outputFile is not sys.stdout:
            self.args.outputFile.close()

    def __align(self, referenceMap: OpticalMap, queryMap: OpticalMap):
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator)
        primaryCorrelationReverse = queryMap.getInitialAlignment(referenceMap, self.primaryGenerator,
                                                                 reverseStrand=True)
        bestPrimaryCorrelation = sorted([primaryCorrelation, primaryCorrelationReverse], key=lambda c: c.getScore())[-1]

        if not bestPrimaryCorrelation.peakPositions.any():
            return None

        secondaryCorrelation = bestPrimaryCorrelation.refine(self.secondaryGenerator, self.args.minAdjustment,
                                                             self.args.peakHeightThreshold)
        if self.args.diagnostics:
            fig = plotCorrelation(bestPrimaryCorrelation)
            fig.savefig(
                f"C:/Users/arcis/Documents/code/optical_maps/output_alignments/correlation_primary_{self.args.queryIds[0]}.svg",
                bbox_inches='tight', pad_inches=0)

            fig = plotRefinedCorrelation(bestPrimaryCorrelation, secondaryCorrelation)
            fig.savefig(
                f"C:/Users/arcis/Documents/code/optical_maps/output_alignments/correlation_refined_{self.args.queryIds[0]}.svg",
                bbox_inches='tight', pad_inches=0)

        alignmentResultRow = self.aligner.align(referenceMap, queryMap, secondaryCorrelation.peakPositions,
                                                secondaryCorrelation.reverseStrand)
        return alignmentResultRow


if __name__ == '__main__':
    main()
