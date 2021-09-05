# %%

from collections import Counter

from matplotlib import cycler  # type: ignore
from matplotlib import rcParams  # type: ignore
from scipy.optimize import minimize_scalar

from src.cmap_reader import CmapReader
from src.optical_map import VectorisedOpticalMap
from src.plot import plotCorrelation
from src.sequence_generator import SequenceGenerator

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def print1To0Ratio(sequence):
    counts = Counter(sequence)
    print(f"1: {counts[1]:,}, 0: {counts[0]:,}")
    print("1 to 0 ratio: %.3f" % ((counts[1]/counts[0])))

# %%


def pipeline(resolution=43, blurRadius=2, plot=False):
    resolution = round(resolution)
    blurRadius = round(blurRadius)
    sequenceGenerator = SequenceGenerator(resolution, blurRadius)
    reader = CmapReader(sequenceGenerator)
    referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
    reference = reader.readReference(referenceFile, 1)

    queryFile = "../data/EXP_REFINEFINAL1.cmap"
    moleculeIds = [171, ]  # [11, 12, 21, 22, 31, 32]

    queries = reader.readQueries(queryFile, moleculeIds)
    query: VectorisedOpticalMap = queries[0]
    query.__reverse()
    result = query.correlate(reference)

    if plot:
        fig = plotCorrelation(result, resolution, False, (51753149, 60405486))
        fig.savefig(f"../plots_irys/plot_molecule{query.moleculeId}_res{resolution}_blur{blurRadius}.svg",
                    bbox_inches='tight', pad_inches=0)

    print(f"res: {resolution}, blur:{blurRadius}, score: {result.peaks.score}")
    return result
    # return result.quality.reverseScore


# %%
minimize_scalar(pipeline, bounds=(25, 60), method='bounded', options={"disp": True, "maxiter": 25, "xatol": 1})
# %%
minimize_scalar(pipeline, bounds=(0, 10), method='bounded', options={"disp": True, "maxiter": 25, "xatol": 1})

# %%
# minimize(pipeline, [40, 5], method="Powell", bounds=((35, 80), (0, 50)), options={'xtol': 0.01})
# %%
for res in [43, 20, 60]:
    result = pipeline(res, plot=True)
    print(result.peaks.score)
    print(len(result.peaks.peaks))

# %%
for blur in [0, 1, 2, 3, 4]:
    result = pipeline(43, blur, True)
    print(result.peaks.score)
    print(len(result.peaks.peaks))

# zobaczyć czy da sie ukraść dynapic programming np z fandom do końcowego alignmentu

# policzyć istotność piku (może podzielić przez 2 największą wartość),
# na podstawie tego dobrać optymalne resolution i blur

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
