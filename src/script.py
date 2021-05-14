# %%

from matplotlib import pyplot as plt
from files.map_reading import readOpticalMap, readReference
from matplotlib import rcParams, cycler
from visuals.plot import plotCorrelation
from collections import Counter

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])
# %load_ext autoreload
# %autoreload 2
# %%

withError = True
resolution = 10
blurRadius = 100
referenceFile = "../data/ecoli_ref.cmap"
referenceMap = readReference(referenceFile, resolution, blurRadius)
# %%
for moleculeId in range(1, 5):
    normalize = True
    sdataMapFile = f"../data/ecoli_ref_{'with' if withError else 'without'}_error_simulated.sdata"

    simulatedMap = readOpticalMap(sdataMapFile, resolution, blurRadius, moleculeId)

    singleVsRef = simulatedMap.correlate(referenceMap, normalize)

    counts = Counter(simulatedMap.sequence)
    print("1 to 0 ratio: %.3f" % ((counts[1]/counts[0])))

    fig = plotCorrelation(singleVsRef, resolution)
    plt.axis('off')
    fig.savefig(f"../plots/plot_molecule{moleculeId}_res{resolution}_blur{blurRadius}_errors_{withError}_normalize_{normalize}.png",
                bbox_inches='tight', pad_inches=0)


# %%
# puścić na danych z fandom
# puścić fandom lub inne, wziąć mapy z dobrym wynikiem
# zobaczyć czy trzeba preprocessing, poprawiać skalowanie
# znaleźć ładny pipeline, może coś z fandom, albo coś gdzie aligner będzie grał dużą rolę
# zobaczyć czy da sie ukraść dynapic programming np z fandom do końcowego alignmentu
# spróbować więcej zmniejszenia rozdzielczości, mniej rozmycia
# %%
