# %%
from collections import Counter

import pandas
from matplotlib import cycler  # type: ignore
from matplotlib import rcParams  # type: ignore
from matplotlib import pyplot as plt
from numpy.lib.function_base import disp

from cmap_reader import AlignmentReader, CmapReader
from plot import plotCorrelation, plotHeatMap
from sequence_generator import SequenceGenerator
from validator import Validator
from tqdm import tqdm
from random import sample
rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


# %%
alignmentsFile = "../data/EXP_REFINEFINAL1.xmap"
referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
queryFile = "../data/EXP_REFINEFINAL1.cmap"


alignmentReader = AlignmentReader()
alignments = alignmentReader.readAlignments(alignmentsFile)

alignmentsCount = 300
resolutions = [32, 48, 64, 256, 1024, 1536]
blurs = [0, 2, 4, 8, 12]
# %%
isoResolutionResults = []
with tqdm(total=alignmentsCount * len(blurs) * len(resolutions)) as progressBar:
    for resolution in resolutions:
        isoBlurResults = []
        for blur in blurs:
            validator = Validator(resolution)
            sequenceGenerator = SequenceGenerator(resolution, blur)
            reader = CmapReader(sequenceGenerator)

            validCount = 0
            for alignment in sample(alignments, alignmentsCount):
                reference = reader.readReference(referenceFile, alignment.chromosome)
                query = reader.readQuery(queryFile, alignment.queryId)
                result = query.correlate(reference, reverse=alignment.reverseStrand)
                isValid = validator.validate(result, alignment)
                if isValid:
                    validCount += 1
                # print(f"confidence: {alignment.confidence}, score: {result.peaks.score}, valid: {isValid}")
                progressBar.update(1)
            isoBlurResults.append(validCount / alignmentsCount)
        isoResolutionResults.append(isoBlurResults)
    # results.append({
    #     'moleculeId': query.moleculeId,
    #     'referenceId': reference.moleculeId,
    #     'alignmentId': alignment.id,
    #     'overlapsReferenceAlignment': isValid,
    #     'referenceAlignmentConfidence': alignment.confidence,
    #     'score': result.peaks.score})

    # if not isValid:
    #     fig = plotCorrelation(result, resolution, False,
    #                         (alignment.refStartPosition,
    #                         alignment.refEndPosition))
    #     fig.savefig(
    #         f"../plots_alignments/plot_molecule{query.moleculeId}_res{resolution}_score{result.peaks.score}_conf{alignment.confidence}.svg", bbox_inches='tight', pad_inches=0)

plotHeatMap(isoResolutionResults, alignmentsCount, resolutions, blurs)
# %%
# 1000 dobrze zmapowanych sekwencji
#  zbadać parametry - heat map, dobrać zakres parametrów tak żeby było widać spadek, do res * blur * 2 < 1000
# wynik: ile % wyników się pokrywa z ref aligner, druga heatmapa z czasami obliczeń
# potem przefiltrować cząsteczki i parametry, tak żeby zawsze mieć 100%, heat mapa parametry -> średnia/mediana score + odchylenie
# kolejny wykres x - jakość,  y - liczba cząsteczek o tej wartości jakości, kernel density
