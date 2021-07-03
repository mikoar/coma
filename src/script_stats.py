# %%

from os import pipe, read
import re
import itertools
from typing import List
from numpy.lib.function_base import disp
import pandas
from matplotlib import pyplot as plt
from matplotlib import rcParams, cycler  # type: ignore
from cmap_reader import AlignmentReader, CmapReader
from optical_map import OpticalMap
from sequence_generator import SequenceGenerator
from plot import plotCorrelation
from collections import Counter
from scipy.optimize import minimize_scalar, minimize

from validator import Validator

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


# %%
resolution = 120
blurRadius = 2
alignmentsFile = "../data/EXP_REFINEFINAL1.xmap"
referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
queryFile = "../data/EXP_REFINEFINAL1.cmap"
validator = Validator(resolution)

sequenceGenerator = SequenceGenerator(resolution, blurRadius)
reader = CmapReader(sequenceGenerator)

alignmentReader = AlignmentReader()
alignments = alignmentReader.readAlignments(alignmentsFile)

for alignment in alignments[:300]:
    reference = reader.readReference(referenceFile, alignment.chromosome)
    query = reader.readQuery(queryFile, alignment.queryId)
    result = query.correlate(reference, reverse=alignment.reverseStrand)
    isValid = validator.validate(result, alignment)

    print(f"confidence: {alignment.confidence}, score: {result.peaks.score}, valid: {isValid}")
    if not isValid:
        fig = plotCorrelation(result, resolution, False,
                              (alignment.refStartPosition,
                               alignment.refEndPosition))
        fig.savefig(f"../plots_alignments/plot_molecule{query.moleculeId}_score{result.peaks.score}_blur{blurRadius}_conf{alignment.confidence}.svg",
                    bbox_inches='tight', pad_inches=0)


# %%
