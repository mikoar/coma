# %%

from files.map_reading import readOpticalMap, readReference
from visuals.plot import plotCorrelation

# %load_ext autoreload
# %autoreload 2
# %%

withError = False
resolution = 10
blurRadius = 100
referenceFile = "../data/ecoli_ref.cmap"
referenceMap = readReference(referenceFile, resolution, blurRadius)
# %%
moleculeId = 2
normalize = True
sdataMapFile = f"../data/ecoli_ref_{'with' if withError else 'without'}_error_simulated.sdata"

simulatedMap = readOpticalMap(sdataMapFile, resolution, blurRadius, moleculeId)

singleVsRef = simulatedMap.correlate(referenceMap, normalize)

fig = plotCorrelation(singleVsRef, resolution)
fig.savefig(f"../plots/plot_molecule{moleculeId}_res{resolution}_blur{blurRadius}_errors_{withError}_normalize_{normalize}.svg",
            bbox_inches='tight', pad_inches=0)


# %%
# TODO: opis pracy: mapowanie map na genom, jeżeli się uda to więcej: składanie odczytów, korekcja błędów, napisać co będzie wynikiem pracy, co zaimplementuję, jak przetestuję, coś co pozwoli komisji ocenić złożóność pracy, nie obiecywać zbyt wiele
# TODO: rozciąganie cząsteczek - do przemyślenia
