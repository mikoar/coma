from src.alignment.aligner import AlignerEngine, Aligner, ChainBuilder
from src.alignment.alignment_position_scorer import AlignmentPositionScorer
from src.alignment.chain_scorer import ChainScorer
from src.alignment.segment_chainer import SegmentChainer, SequentialityScorer
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentConflictResolver
from src.alignment.segments_factory import AlignmentSegmentsFactory
from src.args import Args
from src.correlation.peaks_selector import PeaksSelector
from src.correlation.sequence_generator import SequenceGenerator
from src.extensions.dispatcher import Dispatcher
from src.multi_pass_workflow_coordinator import _MultiPassWorkflowCoordinator
from src.parsers.xmap_reader import XmapReader
from src.workflow_coordinator import _WorkflowCoordinator


class WorkflowCoordinatorFactory:
    def __init__(self, args: Args, dispatcher: Dispatcher, xmapReader: XmapReader):
        self.args = args
        self.dispatcher = dispatcher
        self.xmapReader = xmapReader

    def create(self):
        primaryGenerator = SequenceGenerator(self.args.primaryResolution, self.args.primaryBlur)
        secondaryGenerator = SequenceGenerator(self.args.secondaryResolution, self.args.secondaryBlur)
        scorer = AlignmentPositionScorer(
            self.args.perfectMatchScore,
            self.args.distancePenaltyMultiplier,
            self.args.unmatchedPenalty,
            self.args.endReachingScore)
        if self.args.sequentialityScore==2:
            chainScorer = ChainScorer(self.args.segmentCombinePenalty, self.args.segmentJoinMultiplier, self.args.minScore, self.args.minSubsequentScore)
            aligner = ChainBuilder(scorer, chainScorer, self.args.secondaryResolution*self.args.secondaryBlur)
        else:
            segmentsFactory = AlignmentSegmentsFactory(self.args.minScore, self.args.breakSegmentThreshold)
            alignerEngine = AlignerEngine(self.args.maxPairDistance)
            alignmentSegmentConflictResolver = AlignmentSegmentConflictResolver(
                SegmentChainer(
                    SequentialityScorer(self.args.segmentJoinMultiplier, self.args.sequentialityScore, self.args.segmentCombinePenalty, self.args.secondaryMargin)))
            aligner = Aligner(scorer, segmentsFactory, alignerEngine, alignmentSegmentConflictResolver)
        if self.args.outputMode == "single":
            return _WorkflowCoordinator(
                self.args, primaryGenerator,
                secondaryGenerator,
                aligner,
                self.dispatcher,
                PeaksSelector(self.args.primaryPeaksCount))
        else:
            return _MultiPassWorkflowCoordinator(
                self.args,
                primaryGenerator,
                secondaryGenerator,
                aligner,
                self.dispatcher,
                PeaksSelector(self.args.primaryPeaksCount),
                self.xmapReader)
