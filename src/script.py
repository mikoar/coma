# %%
import matplotlib.patches as patches
from typing import List
import pandas
import itertools
import re
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from vectorise import vectorise
from scipy.signal import correlate
%load_ext autoreload
%autoreload 2

# %%


def read_cmap(file_path, columns=['CMapId', 'ContigLength', 'NumSites', 'SiteID', 'LabelChannel', 'Position']):
    with open(file_path) as file:
        gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
        header_line = list(itertools.islice(gen, 1))[0].strip()
        names = re.split('\s+', header_line)[1:]

    csv = pandas.read_csv(
        file_path, comment="#", delimiter="\t", names=names)
    return csv[columns]


# java -jar OMTools.jar FastaToOM --fastain ../ecoli.fasta --enzyme BspQI --refmapout ../ecoli.cmap
reference_cmap = read_cmap('../data/ecoli_ref.cmap')

# java -jar OMTools.jar OptMapDataGenerator --refmapin ../ecoli_ref.cmap --optmapout ../ecoli_map.cmap
simulated_cmap = read_cmap('../data/ecoli_map.cmap')

# %%
simulated_single_molecule_cmap = simulated_cmap[simulated_cmap.CMapId == 219]

# %%


def positionsToSequence(data):
    return list(vectorise((position for _, position in data.Position.sort_values().iteritems()), 1))


def cross_corelate(seq1: List[int], seq2: List[int]):
    min_len = min(len(seq1), len(seq2))
    corr = correlate(seq1[0:min_len], seq2[0:min_len])
    return corr


def find_shifts(correlation, treshold, query_length):
    shifts = []
    for i, value in enumerate(correlation):
        if value > treshold:
            shifts.append((i - query_length + 1, value))
    return shifts


simulated_seq = positionsToSequence(simulated_cmap)
simulated_single_molecule_seq = positionsToSequence(simulated_single_molecule_cmap)
reference_seq = positionsToSequence(reference_cmap)

single_vs_all = cross_corelate(simulated_seq, simulated_single_molecule_seq)
plt.plot(single_vs_all)

single_vs_all_shifts = find_shifts(single_vs_all, 10, len(simulated_single_molecule_seq))

# %%


# def find_shifts_cumulated(correlation, treshold, query_length):
#     shifts = []
#     for i, value in enumerate(correlation):
#         if value > treshold:
#             print(i)
#             shifts.append((i - query_length + 1, value))
#     return shifts
# %%


single_vs_ref = cross_corelate(simulated_single_molecule_seq, reference_seq)
fig = plt.figure(figsize=(40, 5))
ax = fig.add_axes([0, 0, 1, 1])
ax.plot(single_vs_ref)

# %%
