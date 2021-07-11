# %%
from typing import List
from matplotlib import cycler  # type: ignore
from matplotlib import rcParams
import numpy as np  # type: ignore

from cmap_reader import Alignment, AlignmentReader, CmapReader
from optical_map import OpticalMap

from plot import plotHeatMap
from sequence_generator import SequenceGenerator
from worker import workerFunction
from tqdm import tqdm
from random import sample
import multiprocessing
rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def chunk(list, size):
    for i in range(0, len(list), size):
        yield list[i:i + size]


def getWorkerInputs(alignments: List[Alignment], reference: np.ndarray, queries: List[OpticalMap], resolution: int):
    for alignment in alignments:
        yield (alignment,
               reference,
               next(q for q in queries if q.moleculeId == alignment.queryId),
               resolution)


if __name__ == '__main__':
    multiprocessing.freeze_support()

    alignmentsFile = "../data/EXP_REFINEFINAL1.xmap"
    referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
    queryFile = "../data/EXP_REFINEFINAL1.cmap"

    alignmentsCount = 10
    resolutions = [32, 48, 64, 256, 1024, 1536]
    blurs = [0, 2, 4, 8, 12]

    alignmentReader = AlignmentReader()
    alignments = alignmentReader.readAlignments(alignmentsFile)
    # %%

    isoResolutionResults = []
    with multiprocessing.Pool(processes=6, maxtasksperchild=1) as pool:

        with tqdm(total=alignmentsCount * len(blurs) * len(resolutions)) as progressBar:
            for resolution in resolutions:
                isoBlurResults = []
                for blur in blurs:
                    sequenceGenerator = SequenceGenerator(resolution, blur)
                    reader = CmapReader(sequenceGenerator)

                    validCount = 0
                    sampledAlignments = sample(alignments, alignmentsCount)
                    referenceIds = set(map(lambda a: a.referenceId, sampledAlignments))
                    alignmentsGroupedByReference = [[a for a in sampledAlignments if a.referenceId == r] for r in referenceIds]
                    for alignmentsForReference, referenceId in zip(alignmentsGroupedByReference, referenceIds):
                        reference = reader.readReference(referenceFile, referenceId)
                        queries = reader.readQueries(queryFile, list(map(lambda a: a.queryId, alignmentsForReference)))
                        progressBar.set_description(
                            f"Resolution: {resolution}, blur: {blur}, {len(queries)} queries for reference {referenceId}")

                        poolResults = pool.map(workerFunction, getWorkerInputs(alignmentsForReference, reference.sequence, queries, resolution))
                        validCount += sum(poolResults)
                        progressBar.update(len(alignmentsForReference))

                    isoBlurResults.append(validCount / len(sampledAlignments))

                isoResolutionResults.append(isoBlurResults)

    plotHeatMap(isoResolutionResults, alignmentsCount, blurs, resolutions)


# %%
