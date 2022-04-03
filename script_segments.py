import os
from functools import reduce

from src.alignment.aligner import Aligner
from src.alignment.alignment_results import AlignmentResults
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader

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

    alignmentResultFile = f"output_alignments/alignment_segments_{title}.xmap"
    os.makedirs("output_alignments", exist_ok=True)

    xmapReader = XmapReader()
    cmapReader = CmapReader()

    reference = cmapReader.readReference(referenceFile, referenceId)
    query = cmapReader.readQuery(queryFile, queryId)

    initialGenerator = SequenceGenerator(initialResolution, initialBlur)
    initialCorrelation = query.getInitialAlignment(reference, initialGenerator)
    initialAligner = Aligner(2 * initialResolution * initialBlur)
    nonRefinedAlignmentResult = initialAligner.align(reference, query, initialCorrelation.maxPeak.position)

    refineGenerator = SequenceGenerator(refinedResolution, refineBlur)
    refinedCorrelation = initialCorrelation.refine(refineGenerator)
    print(len(list(refinedCorrelation.peaks)))

    refinedAligner = Aligner(2 * initialResolution * initialBlur)
    refinedAlignmentResults = [refinedAligner.align(reference, query, peak.position) for peak in
                               refinedCorrelation.peaks]
    perfectMatchScore = 400
    scoreMultiplier = 1
    unmatchedPenalty = -100
    minScore = 4 * perfectMatchScore
    breakSegmentThreshold = 600
    filteredAlignmentResults = [
        AlignmentSegments.filterSegments(row, perfectMatchScore, scoreMultiplier, unmatchedPenalty, minScore,
                                         breakSegmentThreshold) for row in refinedAlignmentResults]

    mergedResult = reduce(lambda row1, row2: row1.merge(row2), filteredAlignmentResults)

    refAlignment = xmapReader.readAlignments(alignmentsFile, queryIds=[queryId])[0]
    alignmentResults = AlignmentResults(referenceFile, queryFile,
                                        [
                                            # refAlignment, nonRefinedAlignmentResult,
                                            mergedResult])
    xmapReader.writeAlignments(alignmentResultFile, alignmentResults)

    for result in refinedAlignmentResults + filteredAlignmentResults + [mergedResult]:
        # shifts = list(map(lambda pair: f"{pair.queryShift:.0f}", result.alignedPairs))
        # print(shifts)
        print(result.positions)
