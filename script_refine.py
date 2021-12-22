import os
from statistics import mean

from matplotlib import cycler
from matplotlib import rcParams

from src.alignment.aligner import Aligner
from src.alignment.alignment_results import AlignmentResults
from src.correlation.optical_map import CorrelationResult
from src.correlation.plot import plotCorrelation
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def plot(correlation: CorrelationResult, resolution: int, stage: str, title: str):
    fig = plotCorrelation(correlation, resolution)
    fig.suptitle(f'{stage} correlation {queryId}')
    fig.savefig(f"output_alignments/alignment_refinement_{stage}_{title}.svg", bbox_inches='tight', pad_inches=0)


if __name__ == '__main__':
    referenceId = 21
    queryId = 63655

    baseDir = f'.local_data/NA12878_BSPQI_pipeline_results/output/contigs/alignmolvref/merge/alignmolvref_contig{referenceId}'
    alignmentsFile = f"{baseDir}.xmap"
    referenceFile = f"{baseDir}_r.cmap"
    queryFile = f"{baseDir}_q.cmap"

    initialResolution = 256
    initialBlur = 4
    refinedResolution = 40
    refineBlur = 2
    title = f"{referenceId}_{queryId}_{initialResolution}-{refinedResolution}_{initialBlur}-{refineBlur}"

    alignmentResultFile = f"output_alignments/alignment_refinement_{title}.xmap"
    os.makedirs("output_alignments", exist_ok=True)

    xmapReader = XmapReader()
    cmapReader = CmapReader()

    reference = cmapReader.readReference(referenceFile, referenceId)
    query = cmapReader.readQuery(queryFile, queryId)

    initialGenerator = SequenceGenerator(initialResolution, initialBlur)
    initialCorrelation = query.getInitialAlignment(reference, initialGenerator)
    initialAligner = Aligner(2 * initialResolution * initialBlur)
    nonRefinedAlignmentResult = initialAligner.align(reference, query, initialCorrelation.maxPeak.position)
    plot(initialCorrelation, initialResolution, "initial", title)

    refineGenerator = SequenceGenerator(refinedResolution, refineBlur)
    maxAdjustment = 15000
    refinedCorrelation = initialCorrelation.refine(refineGenerator, maxAdjustment)
    print(len(list(refinedCorrelation.peaks)))

    refinedAligner = Aligner(2 * initialResolution * initialBlur)
    plot(refinedCorrelation, refinedResolution, "refined", title)
    refinedAlignmentResults = [refinedAligner.align(reference, query, peak.position) for peak in
                               refinedCorrelation.peaks]

    refAlignment = xmapReader.readAlignments(alignmentsFile, queryIds=[queryId])[0]
    alignmentResults = AlignmentResults(referenceFile, queryFile,
                                        [refAlignment, nonRefinedAlignmentResult] + refinedAlignmentResults)
    xmapReader.writeAlignments(alignmentResultFile, alignmentResults)

    for result in refinedAlignmentResults:
        shifts = list(map(lambda pair: pair.queryShift, result.alignedPairs))
        print(shifts)
        print(mean(shifts))
