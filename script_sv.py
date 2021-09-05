# %%

from os import pipe
import re
import itertools
from typing import List
from numpy.lib.function_base import disp
import pandas
from matplotlib import pyplot as plt
from matplotlib import rcParams, cycler  # type: ignore
from cmap_reader import CmapReader
from optical_map import VectorisedOpticalMap
from sequence_generator import SequenceGenerator
from plot import plotCorrelation
from collections import Counter
from scipy.optimize import minimize_scalar, minimize

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])

testCases = [{"moleculeId": 111, "reverse": True, "referenceId": 1, "referenceRanges":
              [(2694767, 13467866), (13689849, 16979795), (16892657, 16979795), (16892657, 16979795),
               (16960960, 17057644), (17080656, 17122165), (17185676, 83904286), (78482505, 78496448)]},
             {"moleculeId": 2252, "reverse": False, "referenceId": 24, "referenceRanges":
              [(9861441, 10029352), (13243421, 13420247)]},
             {"moleculeId": 2060, "reverse": True, "referenceId": 23, "referenceRanges":
              [(51606077, 51956334), (51663121, 51814969)]}]  # one is reverse, one is not


def pipeline(resolution=43, blurRadius=2, plot=False):
    testCase = testCases[2]
    dataPath = "../.local_data/GM09888_pipeline_results.tar/output/contigs/exp_refineFinal1_sv"
    referenceFile = dataPath + "/EXP_REFINEFINAL1_r.cmap"
    resolution = round(resolution)
    blurRadius = round(blurRadius)
    sequenceGenerator = SequenceGenerator(resolution, blurRadius)
    reader = CmapReader(sequenceGenerator)
    reference = reader.readReference(referenceFile, testCase["referenceId"])

    queryFile = dataPath + "/EXP_REFINEFINAL1_q.cmap"
    moleculeIds = [testCase["moleculeId"]]

    queries = reader.readQueries(queryFile, moleculeIds)
    query: VectorisedOpticalMap = queries[0]
    if testCase["reverse"]:
        query.__reverse()
    result = query.correlate(reference)

    if plot:
        referenceRanges = testCase["referenceRanges"]
        fig = plotCorrelation(result, resolution, False, referenceRanges)
        fig.savefig(f"../plots_irys/plot_molecule{query.moleculeId}_res{resolution}_blur{blurRadius}.svg",
                    bbox_inches='tight', pad_inches=0)

    print(f"res: {resolution}, blur:{blurRadius}, score: {result.peaks.score}")
    return result
    # return result.quality.reverseScore


pipeline(plot=True)

# %%
