from __future__ import annotations

import argparse
import sys
from typing import TextIO, List, NamedTuple

from p_tqdm import p_imap

from src.alignment.aligner import Aligner, AlignerEngine
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResults
from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.correlation.optical_map import OpticalMap
from src.correlation.plot import plotRefinedCorrelation
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader


def main():
    parser = argparse.ArgumentParser(description="Optical map aligner.")
    parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"))
    parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"))
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("-r1", "--primaryResolution", dest="primaryResolution", type=int, default=256)
    parser.add_argument("-b1", "--primaryBlur", dest="primaryBlur", type=int, default=4)
    parser.add_argument("-r2", "--secondaryResolution", dest="secondaryResolution", type=int, default=40)
    parser.add_argument("-b2", "--secondaryBlur", dest="secondaryBlur", type=int, default=2)
    parser.add_argument("-ma", "--minAdjustment", dest="minAdjustment", type=int, default=1024)
    parser.add_argument("-rId", "--referenceIDs", dest="referenceIds", type=int, nargs="*")
    parser.add_argument("-qId", "--queryIDs", dest="queryIds", type=int, nargs="*")
    parser.add_argument("-c", "--cpus", dest="numberOfCpus", type=int, default=None)
    parser.add_argument("-d", "--maxDistance", dest="maxDistance", type=int, default=1000)
    parser.add_argument("-sp", "--perfectMatchScore", dest="perfectMatchScore", type=int, default=800)
    parser.add_argument("-sm", "--scoreMultiplier", dest="scoreMultiplier", type=float, default=1.)
    parser.add_argument("-su", "--unmatchedPenalty", dest="unmatchedPenalty", type=int, default=-100)
    parser.add_argument("-ms", "--minScore", dest="minScore", type=int, default=1600)
    parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=600)
    parser.add_argument("-p", "--plot", dest="plot", action="store_true")

    args: Args = parser.parse_args()  # type: ignore
    Program(args).run()


class Args(NamedTuple):
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO
    primaryResolution: int
    primaryBlur: int
    secondaryResolution: int
    secondaryBlur: int
    minAdjustment: int
    referenceIds: List[int]
    queryIds: List[int]
    numberOfCpus: int | None
    maxDistance: int
    perfectMatchScore: int
    scoreMultiplier: int
    unmatchedPenalty: int
    minScore: int
    breakSegmentThreshold: int
    plot: bool


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

        secondaryCorrelation = bestPrimaryCorrelation.refine(self.secondaryGenerator, self.args.minAdjustment)
        if self.args.plot:
            fig = plotRefinedCorrelation(bestPrimaryCorrelation, secondaryCorrelation, self.args.secondaryResolution)
            fig.savefig(
                f"C:/Users/arcis/Documents/code/optical_maps/output_alignments/correlation_refined_{self.args.queryIds[0]}.svg",
                bbox_inches='tight', pad_inches=0)

        alignmentResultRow = self.aligner.align(referenceMap, queryMap, secondaryCorrelation.peakPositions,
                                                secondaryCorrelation.reverseStrand)
        return alignmentResultRow


if __name__ == '__main__':
    main()
