import itertools

import numpy as np
from pandas import DataFrame

from diagnostic.plot import plotHeatMap

if __name__ == '__main__':
    fileDirectory = "../.local_data/simulated/"
    queryName = f"filtered_molecules_1500_200T_scale001_meas500_res600_fpr000002"
    referenceFile = fileDirectory + "chr1_r_unchanged.cmap"
    resolutions = [50, 100, 200, 400, 800]
    blurs = [1, 2, 3, 4, 5]

    statisticsPerParam = []
    for resolution in resolutions:
        for blur in blurs:
            alignmentCompareFile = fileDirectory + f"output/2cc_params/{queryName}_r2_{resolution}_b2_{blur}.tsv"

            with open(alignmentCompareFile, "r") as f:
                headers = [li.replace("# ", "").split("\t") for li in itertools.islice(f, 7)]

            statistics = {h[0]: float(h[1]) for h in headers}

            total = statistics["Overlapping"] + statistics["NonOverlapping"] + statistics["FirstOnly"] + statistics["SecondOnly"]
            overlappingRate = statistics["Overlapping"] / total
            identity = statistics["AvgOverlappingIdentity"] * overlappingRate
            statisticsPerParam.append({"Resolution": resolution, "Blur": blur, "Identity": identity})

    df = DataFrame(statisticsPerParam)
    plotHeatMap(
        np.array(df["Identity"]).reshape(len(blurs), len(resolutions)),
        f"../.local_data/simulated/output/{queryName}_secondary_heatmap.svg",
        blurs,
        resolutions,
        title=f"")
