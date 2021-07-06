# %%
from collections import Counter

import pandas
from matplotlib import cycler  # type: ignore
from matplotlib import rcParams  # type: ignore
from matplotlib import pyplot as plt
from numpy.lib.function_base import disp

from cmap_reader import AlignmentReader, CmapReader
from plot import plotCorrelation
from sequence_generator import SequenceGenerator
from validator import Validator

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


# %%
resolution = 750
blurRadius = 2
alignmentsFile = "../data/EXP_REFINEFINAL1.xmap"
referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
queryFile = "../data/EXP_REFINEFINAL1.cmap"


validator = Validator(resolution)

sequenceGenerator = SequenceGenerator(resolution, blurRadius)
reader = CmapReader(sequenceGenerator)

alignmentReader = AlignmentReader()
alignments = alignmentReader.readAlignments(alignmentsFile)

results = []

for alignment in alignments[:300]:
    reference = reader.readReference(referenceFile, alignment.chromosome)
    query = reader.readQuery(queryFile, alignment.queryId)
    result = query.correlate(reference, reverse=alignment.reverseStrand)
    isValid = validator.validate(result, alignment)

    print(f"confidence: {alignment.confidence}, score: {result.peaks.score}, valid: {isValid}")

    results.append({
        'moleculeId': query.moleculeId,
        'referenceId': reference.moleculeId,
        'alignmentId': alignment.id,
        'overlapsReferenceAlignment': isValid,
        'referenceAlignmentConfidence': alignment.confidence,
        'score': result.peaks.score})

    if not isValid:
        fig = plotCorrelation(result, resolution, False,
                              (alignment.refStartPosition,
                               alignment.refEndPosition))
        fig.savefig(
            f"../plots_alignments/plot_molecule{query.moleculeId}_res{resolution}_score{result.peaks.score}_conf{alignment.confidence}.svg", bbox_inches='tight', pad_inches=0)


# %%
# 1000 dobrze zmapowanych sekwencji
#  zbadać parametry - heat map, dobrać zakres parametrów tak żeby było widać spadek, do res * blur * 2 < 1000
# wynik: ile % wyników się pokrywa z ref aligner, druga heatmapa z czasami obliczeń
# potem przefiltrować cząsteczki i parametry, tak żeby zawsze mieć 100%, heat mapa parametry -> średnia/mediana score + odchylenie
# kolejny wykres x - jakość,  y - liczba cząsteczek o tej wartości jakości, kernel density
