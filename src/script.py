# %%

import re
import itertools
import pandas
from matplotlib import pyplot as plt
from matplotlib import rcParams, cycler
from files.cmap_reader import CmapReader
from maps.optical_map import OpticalMap
from maps.sequence_generator import SequenceGenerator
from visuals.plot import plotCorrelation
from collections import Counter

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])
# %%

resolution = 10
blurRadius = 30
normalize = True
sequenceGenerator = SequenceGenerator(resolution, blurRadius)
reader = CmapReader(sequenceGenerator)
referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
queryFile = "../data/EXP_REFINEFINAL1.cmap"
moleculeIds = [12451, ]  # [11, 12, 21, 22, 31, 32]

reference = reader.readReference(referenceFile)
queries = reader.readQueryMaps(queryFile, moleculeIds)

query: OpticalMap
for query in queries:
    # query.reverse()
    result = query.correlate(reference, normalize)
    fig = plotCorrelation(result, resolution)
    fig.savefig(f"../plots_irys/plot_molecule{query.moleculeId}_res{resolution}_blur{blurRadius}_normalize_{normalize}.svg",
                bbox_inches='tight', pad_inches=0)

    counts = Counter(query.sequence)
    print("1 to 0 ratio: %.3f" % ((counts[1]/counts[0])))


# %%
# puścić na danych z fandom
# puścić fandom lub inne, wziąć mapy z dobrym wynikiem
# zobaczyć czy trzeba preprocessing, poprawiać skalowanie
# znaleźć ładny pipeline, może coś z fandom, albo coś gdzie aligner będzie grał dużą rolę
# zobaczyć czy da sie ukraść dynapic programming np z fandom do końcowego alignmentu
# spróbować więcej zmniejszenia rozdzielczości, mniej rozmycia
# %%


# def __getCmapColumnNames(filePath):
#     with open(filePath) as file:
#         gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
#         header_line = list(itertools.islice(gen, 1))[0].strip()
#         names = re.split('\s+', header_line)[1:]
#     return names


# filePath = "../data/EXP_REFINEFINAL1.cmap"
# maps = pandas.read_csv(
#     filePath, comment="#", delimiter="\t", names=__getCmapColumnNames(filePath))

# # %%
