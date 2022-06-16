import argparse
import sys
from typing import TextIO, List

from tqdm import tqdm

from src.alignment.aligner import Aligner, AlignerEngine
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.alignment_results import AlignmentResults
from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.correlation.optical_map import OpticalMap
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader

if __name__ == '__main__':
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
    parser.add_argument("-c", "--cpus", dest="numberOfCpus", type=int, default=1)
    parser.add_argument("-d", "--maxDistance", dest="maxDistance", type=int, default=2048)
    parser.add_argument("-sp", "--perfectMatchScore", dest="perfectMatchScore", type=int, default=400)
    parser.add_argument("-sm", "--scoreMultiplier", dest="scoreMultiplier", type=float, default=1.)
    parser.add_argument("-su", "--unmatchedPenalty", dest="unmatchedPenalty", type=int, default=-100)
    parser.add_argument("-ms", "--minScore", dest="minScore", type=int, default=1600)
    parser.add_argument("-bs", "--breakSegmentThreshold", dest="breakSegmentThreshold", type=int, default=600)

    args = parser.parse_args()
    referenceFile: TextIO = args.referenceFile
    queryFile: TextIO = args.queryFile
    outputFile: TextIO = args.outputFile
    primaryResolution: int = args.primaryResolution
    primaryBlur: int = args.primaryBlur
    secondaryResolution: int = args.secondaryResolution
    secondaryBlur: int = args.secondaryBlur
    minAdjustment: int = args.minAdjustment
    referenceIds: List[int] = args.referenceIds
    queryIds: List[int] = args.queryIds
    numberOfCpus: int = args.numberOfCpus
    maxDistance: int = args.maxDistance
    perfectMatchScore: int = args.perfectMatchScore
    scoreMultiplier: int = args.scoreMultiplier
    unmatchedPenalty: int = args.unmatchedPenalty
    minScore: int = args.minScore
    breakSegmentThreshold: int = args.breakSegmentThreshold

    cmapReader = CmapReader()
    xmapReader = XmapReader()
    primaryGenerator = SequenceGenerator(primaryResolution, primaryBlur)
    secondaryGenerator = SequenceGenerator(secondaryResolution, secondaryBlur)
    scorer = AlignmentPositionScorer(perfectMatchScore, scoreMultiplier, unmatchedPenalty)
    segmentsFactory = AlignmentSegmentsFactory(minScore, breakSegmentThreshold)
    alignerEngine = AlignerEngine(maxDistance)
    alignmentSegmentConflictResolver = AlignmentSegmentConflictResolver(SegmentChainer())
    aligner = Aligner(scorer, segmentsFactory, alignerEngine, alignmentSegmentConflictResolver)

    referenceMaps: List[OpticalMap]
    queryMaps: List[OpticalMap]
    with referenceFile:
        referenceMaps = cmapReader.readReferences(referenceFile, referenceIds)
    with queryFile:
        queryMaps = cmapReader.readQueries(queryFile, queryIds)


    def align(referenceMap: OpticalMap, queryMap: OpticalMap):
        primaryCorrelation = queryMap.getInitialAlignment(referenceMap, primaryGenerator)
        primaryCorrelationReverse = queryMap.getInitialAlignment(referenceMap, primaryGenerator, reverseStrand=True)
        bestPrimaryCorrelation, isReverse = \
            sorted([(primaryCorrelation, False), (primaryCorrelationReverse, True)], key=lambda c: c[0].getScore())[-1]

        if not bestPrimaryCorrelation.peakPositions.any():
            return None

        secondaryCorrelation = bestPrimaryCorrelation.refine(secondaryGenerator, minAdjustment)
        alignmentResultRow = aligner.align(referenceMap, queryMap, secondaryCorrelation.peakPositions, isReverse)
        return alignmentResultRow


    count = len(referenceMaps) * len(queryMaps)
    alignmentResultRows = [r for r in tqdm((align(r, q) for r in referenceMaps for q in queryMaps), total=count) if
                           r is not None]
    alignmentResult = AlignmentResults(referenceFile.name, queryFile.name, alignmentResultRows)

    xmapReader.writeAlignments(outputFile, alignmentResult)
    if outputFile is not sys.stdout:
        outputFile.close()
