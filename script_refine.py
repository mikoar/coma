import os
from functools import reduce

from matplotlib import cycler
from matplotlib import rcParams
from matplotlib.figure import Figure

from src.alignment.aligner import Aligner
from src.alignment.alignment_results import AlignmentResults
from src.correlation.plot import plotRefinedCorrelation, plotCorrelation
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def savePlot(fig: Figure, suptitle: str, title: str):
    fig.suptitle(f'{suptitle} correlation, chromosome {referenceId}, molecule {queryId}')
    fig.savefig(f"output_alignments/alignment_refinement_{suptitle}_{title}.svg", bbox_inches='tight', pad_inches=0)


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
    savePlot(plotCorrelation(initialCorrelation, initialResolution), "initial", title)

    refineGenerator = SequenceGenerator(refinedResolution, refineBlur)
    refinedCorrelation = initialCorrelation.refine(refineGenerator)
    print(len(list(refinedCorrelation.peaks)))

    refinedAligner = Aligner(2 * initialResolution * initialBlur)
    savePlot(plotRefinedCorrelation(initialCorrelation, refinedCorrelation, refinedResolution), "refined", title)
    refinedAlignmentResults = [refinedAligner.align(reference, query, peak.position) for peak in
                               refinedCorrelation.peaks]
    mergedResult = reduce(lambda row1, row2: row1.merge(row2), refinedAlignmentResults)

    refAlignment = xmapReader.readAlignments(alignmentsFile, queryIds=[queryId])[0]
    alignmentResults = AlignmentResults(referenceFile, queryFile,
                                        [refAlignment, nonRefinedAlignmentResult, mergedResult])
    xmapReader.writeAlignments(alignmentResultFile, alignmentResults)

    for result in refinedAlignmentResults + [mergedResult]:
        shifts = list(map(lambda pair: f"{pair.queryShift:.0f}", result.alignedPairs))
        print(shifts)
        print(result.alignedPairs)
