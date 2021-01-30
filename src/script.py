# %%
from typing import List
import pandas
import itertools
import re
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
from detect import detect
from scipy.signal import correlate
# %%

file_path = '../../data/mLemCat1_Saphyr_DLE1_2857687.cmap'


def read_cmap(file_path):
    with open(file_path) as file:
        gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
        header_line = list(itertools.islice(gen, 1))[0].strip()
        names = re.split('\s+', header_line)[1:]

    csv = pandas.read_csv(
        file_path, comment="#", delimiter="\t", names=names)
    return csv


csv = read_cmap(file_path)
# %%
map1 = csv[csv.CMapId == 29]
map2 = csv[csv.CMapId == 31]

# %%


def positionsToSequence(data):
    return list(detect((position for _, position in data.Position.iteritems())))


def cross_corelate(seq1: List[int], seq2: List[int]):
    min_len = min(len(seq1), len(seq2))
    corr = correlate(seq1[0:min_len], seq2[0:min_len])
    return corr


seq1 = positionsToSequence(map1)
seq2 = positionsToSequence(map2)

corr = cross_corelate(seq1, seq2)
plt.plot(corr)

# %%

corr2 = cross_corelate(seq1, seq1)
plt.plot(corr2)
# %%
