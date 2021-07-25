# %%
from random import Random
from typing import Dict, List

import numpy as np
import pandas as pd
from matplotlib import cycler  # type: ignore
from matplotlib import rcParams
from p_tqdm import p_map
from tqdm import tqdm

from src.alignment import Alignment
from src.cmap_reader import AlignmentReader, CmapReader
from src.optical_map import OpticalMap, Peaks
from src.plot import plotHeatMap
from src.sequence_generator import SequenceGenerator
from src.validator import Validator

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def getWorkerInputs(alignments: List[Alignment], reference: np.ndarray, queries: List[OpticalMap], resolution: int):
    for alignment in alignments:
        yield (alignment,
               reference,
               next(q for q in queries if q.moleculeId == alignment.queryId),
               resolution)


def alignmentsToDict(a: Alignment, score: float, resolution: int, blur: int, isValid: bool):
    return {
        'resolution': resolution,
        'blur': blur,
        'alignmentId': a.id,
        'queryId': a.queryId,
        'referenceId': a.referenceId,
        'confidence': a.confidence,
        'score': score,
        'isValid': isValid
    }


indexCols = ['resolution', 'blur', 'alignmentId']


def initAlignmentsFile(file):
    pd.DataFrame(columns=[
        'resolution',
        'blur',
        'alignmentId',
        'queryId',
        'referenceId',
        'confidence',
        'score',
        'isValid'
    ]).set_index(indexCols).to_csv(file, mode='w')


def appendAlignmentsToFile(alignments: List[Dict], file):
    pd.DataFrame(alignments).set_index(indexCols).to_csv(file, mode='a', header=False)


def alignWithReference(input):
    (alignment, reference, query, resolution) = input
    result = query.correlate(reference, reverseStrand=alignment.reverseStrand)
    validator = Validator(resolution)
    peaks = Peaks(result)
    isMaxPeakValid = validator.validate(peaks.max, alignment)

    return (1 if isMaxPeakValid else 0, peaks.getRelativeScore(alignment, validator))
# %%


if __name__ == '__main__':
    alignmentsFile = "data/EXP_REFINEFINAL1.xmap"
    referenceFile = "data/hg19_NT.BSPQI_0kb_0labels.cmap"
    queryFile = "data/EXP_REFINEFINAL1.cmap"

    df = pd.DataFrame()

    alignmentsCount = 1000
    resolutions = [32, 48, 64, 128, 256, 512]
    blurs = [0, 2, 4, 8, 16]
    title = f"count_{alignmentsCount}_res_{','.join(str(x) for x in resolutions)}_blur_{','.join(str(x) for x in blurs)}"

    alignmentsResultFile = f"output_heatmap/result_{title}.csv"
    initAlignmentsFile(alignmentsResultFile)

    alignmentReader = AlignmentReader()
    alignments = alignmentReader.readAlignments(alignmentsFile)
    # %%

    isoResolutionResults = []
    with tqdm(total=alignmentsCount * len(blurs) * len(resolutions)) as progressBar:
        for resolution in resolutions:
            isoBlurResults = []
            for blur in blurs:
                sequenceGenerator = SequenceGenerator(resolution, blur)
                reader = CmapReader(sequenceGenerator)

                validAlignments = []

                validCount = 0
                sampledAlignments = Random(123).sample(alignments, alignmentsCount,)
                referenceIds = set(map(lambda a: a.referenceId, sampledAlignments))
                alignmentsGroupedByReference = [[a for a in sampledAlignments if a.referenceId == r] for r in referenceIds]
                for alignmentsForReference, referenceId in zip(alignmentsGroupedByReference, referenceIds):
                    reference = reader.readReference(referenceFile, referenceId)
                    queries = reader.readQueries(queryFile, list(map(lambda a: a.queryId, alignmentsForReference)))
                    progressBar.set_description(
                        f"Resolution: {resolution}, blur: {blur}, {len(queries)} queries for reference {referenceId}")

                    poolResults = p_map(alignWithReference, list(getWorkerInputs(alignmentsForReference, reference.sequence, queries, resolution)), num_cpus=8)
                    areValid, scores = zip(*poolResults)

                    alignmentDataToStore = [alignmentsToDict(a, score, resolution, blur, isValid)
                                            for a, score, isValid in zip(alignmentsForReference, scores, areValid)]

                    appendAlignmentsToFile(alignmentDataToStore, alignmentsResultFile)

                    validCount += sum(areValid)
                    progressBar.update(len(alignmentsForReference))

                isoBlurResults.append(validCount / len(sampledAlignments) * 100)

            isoResolutionResults.append(isoBlurResults)

    plotHeatMap(isoResolutionResults, f"output_heatmap/heatmap_{title}.svg", blurs, resolutions)

# %%
# 1000 dobrze zmapowanych sekwencji
#  zbadać parametry - heat map, dobrać zakres parametrów tak żeby było widać spadek, do res * blur * 2 < 1000
# wynik: ile % wyników się pokrywa z ref aligner, druga heatmapa z czasami obliczeń
# potem przefiltrować cząsteczki i parametry, tak żeby zawsze mieć 100%, heat mapa parametry -> średnia/mediana score + odchylenie
# kolejny wykres x - jakość,  y - liczba cząsteczek o tej wartości jakości, kernel density

# ustalić z których danych korzystam


# powtórzyć heat mapę, peak z obszaru oczekiwanego (lub max jeśli brak peaka) porównać z 5. peakiem, wyłączyć min-height
