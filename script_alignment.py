import os
from random import Random
from typing import List

import numpy as np
import pandas as pd
from matplotlib import cycler
from matplotlib import rcParams
from p_tqdm import p_map
from tqdm import tqdm

from src.alignment.aligner import Aligner
from src.alignment.alignment_comparer import AlignmentComparer, AlignmentComparisonResult
from src.alignment.alignment_results import AlignmentResults
from src.correlation.bionano_alignment import BionanoAlignment
from src.correlation.optical_map import OpticalMap
from src.correlation.sequence_generator import SequenceGenerator
from src.correlation.validator import Validator
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def getWorkerInputs(alignments: List[BionanoAlignment], reference: np.ndarray, queries: List[OpticalMap],
                    resolution: int, generator: SequenceGenerator):
    for alignment in alignments:
        yield (alignment,
               reference,
               next(q for q in queries if q.moleculeId == alignment.queryId),
               resolution,
               generator)


def initFile(file):
    pd.DataFrame(columns=[
        "query1Coverage", "query2Coverage", "queryId", "referenceId", "pairs1", "pairs2"
    ]).to_csv(file, mode='w')


def appendToFile(items: List[AlignmentComparisonResult], file):
    pd.DataFrame.from_dict([{
        "query1Coverage": i.query1Coverage,
        "query2Coverage": i.query2Coverage,
        "queryId": i.queryId,
        "referenceId": i.referenceId,
        "pairs1": i.pairs1,
        "pairs2": i.pairs2
    } for i in items if i]).to_csv(file, mode='a', header=False)


def alignWithReference(params):
    refAlignment: BionanoAlignment
    reference: OpticalMap
    query: OpticalMap
    resolution: int
    generator: SequenceGenerator
    refAlignment, reference, query, resolution, generator = params
    peaks = query.getInitialAlignment(reference, generator, refAlignment.reverseStrand)
    validator = Validator(resolution)
    isMaxPeakValid = validator.validate(peaks.maxPeak, refAlignment)
    if not isMaxPeakValid:
        return

    alignmentResultRow = Aligner(3000).align(reference, query, peaks.maxPeak.positionInReference,
                                             refAlignment.reverseStrand)
    return AlignmentComparer().compare(refAlignment, alignmentResultRow), alignmentResultRow


if __name__ == '__main__':
    chromosome = 21
    baseDir = f'.local_data/NA12878_BSPQI_pipeline_results/output/contigs/alignmolvref/merge/alignmolvref_contig{chromosome}'
    alignmentsFile = f"{baseDir}.xmap"
    referenceFile = f"{baseDir}_r.cmap"
    queryFile = f"{baseDir}_q.cmap"
    # alignmentsFile = "data/NA12878_BSPQI/EXP_REFINEFINAL1.xmap"
    # referenceFile = "data/NA12878_BSPQI/hg19_NT.BSPQI_0kb_0labels.cmap"
    # queryFile = "data/NA12878_BSPQI/EXP_REFINEFINAL1.cmap"

    xmapReader = XmapReader()
    alignments = xmapReader.readAlignments(alignmentsFile)
    alignmentsCount = 100  # len(alignments)
    resolution = 256  # [128, 256, 512, 1024]
    blur = 4  # [0, 2, 4, 8, 16]
    title = f"count_{alignmentsCount}_res_{resolution}_blur_{blur}"

    alignmentComparisonResultFile = f"output_alignments/compare_{chromosome}_{title}.csv"
    alignmentResultFile = f"output_alignments/alignment_v2_{chromosome}_{title}.xmap"
    os.makedirs("output_alignments", exist_ok=True)
    initFile(alignmentComparisonResultFile)

    with tqdm(total=alignmentsCount) as progressBar:
        reader = CmapReader()
        alignmentResultRows = []
        sampledAlignments = [a for a in Random(123).sample([a for a in alignments], alignmentsCount)]
        referenceIds = set(map(lambda a: a.referenceId, sampledAlignments))
        alignmentsGroupedByReference = [[a for a in sampledAlignments if a.referenceId == r] for r in
                                        referenceIds]
        for alignmentsForReference, referenceId in zip(alignmentsGroupedByReference, referenceIds):
            reference = reader.readReference(referenceFile, referenceId)
            queries = reader.readQueries(queryFile, list(map(lambda a: a.queryId, alignmentsForReference)))
            progressBar.set_description(
                f"Resolution: {resolution}, blur: {blur}, {len(queries)} queries for reference {referenceId}")

            poolResults = p_map(alignWithReference, list(
                getWorkerInputs(alignmentsForReference, reference, queries, resolution,
                                SequenceGenerator(resolution, blur))), num_cpus=8)
            alignmentComparisonResults, alignmentResultRowsForReference = zip(*[r for r in poolResults if r])
            alignmentResultRows += alignmentResultRowsForReference
            appendToFile(alignmentComparisonResults, alignmentComparisonResultFile)
            progressBar.update(len(alignmentsForReference))

        alignmentResults = AlignmentResults(referenceFile, queryFile, [r for r in alignmentResultRows])
        xmapReader.writeAlignments(alignmentResultFile, alignmentResults)

        # TODO: wyzualizacje porównań alignmentów, zaaplikować ucinanie segmentów, zaznaczyć na wizualizacji gdzie ucina
