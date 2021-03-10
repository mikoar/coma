# %%

from files.map_reading import readOpticalMap, readReference
from visuals.plot import plotCorrelation

# %load_ext autoreload
# %autoreload 2


# %%
withError = False
referenceFile = rf"D:\Dysk Google\bioinf\mgr\data\simulated_{'with' if withError else 'without'}_error\ecoli.cmap"
mapFile = rf"D:\Dysk Google\bioinf\mgr\data\simulated_{'with' if withError else 'without'}_error\ecoli_simulated.cmap"
sdataMapFile = rf"D:\Dysk Google\bioinf\mgr\data\simulated_{'with' if withError else 'without'}_error\ecoli_simulated.sdata"
moleculeId = 10
resolution = 10
blurRadius = 100

simulatedMap = readOpticalMap(sdataMapFile, resolution, blurRadius, moleculeId)
referenceMap = readReference(referenceFile, resolution, blurRadius)

singleVsRef = simulatedMap.correlate(referenceMap)

fig = plotCorrelation(singleVsRef, resolution)
fig.savefig(f"../plots/plot_molecule{moleculeId}_res{resolution}_blur{blurRadius}_errors_{withError}.svg",  bbox_inches='tight', pad_inches=0)
# print(find_shifts(single_vs_ref, max(single_vs_ref)-3, len(simulatedMap)))

# %%
# TODO: znormalizować wykres - podzielić przez wartość oczekiwaną - korelacja z wektorem jedynek o długości odczytu
# TODO: wyzerować wszystkie błędy
# TODO: opis pracy: mapowanie map na genom, jeżeli się uda to więcej: składanie odczytów, korekcja błędów, napisać co będzie wynikiem pracy, co zaimplementuję, jak przetestuję, coś co pozwoli komisji ocenić złożóność pracy, nie obiecywać zbyt wiele
# TODO: rozciąganie cząsteczek - do przemyślenia
