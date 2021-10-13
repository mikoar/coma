import os
from random import Random
from typing import List, Tuple

import numpy as np
import pandas as pd
from matplotlib import cycler
from matplotlib import rcParams
from p_tqdm import p_map
from tqdm import tqdm

from src.alignment.aligner import Aligner
from src.alignment.alignment_comparer import AlignmentComparer, AlignmentComparisonResult
from src.correlation.alignment import Alignment
from src.correlation.cmap_reader import AlignmentReader, CmapReader
from src.correlation.optical_map import VectorisedOpticalMap, Peaks
from src.correlation.sequence_generator import SequenceGenerator
from src.correlation.validator import Validator

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def getWorkerInputs(alignments: List[Alignment], reference: np.ndarray, queries: List[VectorisedOpticalMap],
                    resolution: int):
    for alignment in alignments:
        yield (alignment,
               reference,
               next(q for q in queries if q.moleculeId == alignment.queryId),
               resolution)


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


def alignWithReference(params: Tuple[Alignment, VectorisedOpticalMap, VectorisedOpticalMap, int]):
    refAlignment, reference, query, resolution = params
    result = query.correlate(reference.sequence, reverseStrand=refAlignment.reverseStrand)
    validator = Validator(resolution)
    peaks = Peaks(result)
    isMaxPeakValid = validator.validate(peaks.max, refAlignment)
    if not isMaxPeakValid:
        return

    alignmentResult = Aligner(1000).align(reference, query, peaks.max.positionInReference, refAlignment.reverseStrand)
    return AlignmentComparer().compare(refAlignment, alignmentResult)


if __name__ == '__main__':
    baseDir = '.local_data/NA12878_BSPQI_pipeline_results/output/contigs/alignmolvref/merge/alignmolvref_contig21'
    alignmentsFile = f"{baseDir}.xmap"
    referenceFile = f"{baseDir}_r.cmap"
    queryFile = f"{baseDir}_q.cmap"
    # alignmentsFile = "data/NA12878_BSPQI/EXP_REFINEFINAL1.xmap"
    # referenceFile = "data/NA12878_BSPQI/hg19_NT.BSPQI_0kb_0labels.cmap"
    # queryFile = "data/NA12878_BSPQI/EXP_REFINEFINAL1.cmap"

    alignmentReader = AlignmentReader()
    alignments = alignmentReader.readAlignments(alignmentsFile)
    alignmentsCount = 100  # len(alignments)
    resolution = 256  # [128, 256, 512, 1024]
    blur = 4  # [0, 2, 4, 8, 16]
    title = f"count_{alignmentsCount}_res_{resolution}_blur_{blur}"

    alignmentComparisonResultFile = f"output_alignments/result_{title}.csv"
    os.makedirs("output_alignments", exist_ok=True)
    initFile(alignmentComparisonResultFile)
    # %%

    with tqdm(total=alignmentsCount) as progressBar:
        sequenceGenerator = SequenceGenerator(resolution, blur)
        reader = CmapReader(sequenceGenerator)

        sampledAlignments = [a for a in Random(123).sample([a for a in alignments], alignmentsCount)]
        referenceIds = set(map(lambda a: a.referenceId, sampledAlignments))
        alignmentsGroupedByReference = [[a for a in sampledAlignments if a.referenceId == r] for r in
                                        referenceIds]
        for alignmentsForReference, referenceId in zip(alignmentsGroupedByReference, referenceIds):
            reference = reader.readReference(referenceFile, referenceId)
            queries = reader.readQueries(queryFile, list(map(lambda a: a.queryId, alignmentsForReference)))
            progressBar.set_description(
                f"Resolution: {resolution}, blur: {blur}, {len(queries)} queries for reference {referenceId}")

            alignmentComparisonResults = p_map(alignWithReference, list(
                getWorkerInputs(alignmentsForReference, reference, queries, resolution)), num_cpus=8)

            appendToFile(alignmentComparisonResults, alignmentComparisonResultFile)
            progressBar.update(len(alignmentsForReference))