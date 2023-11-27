import matplotlib
import pandas
from pandas import DataFrame

from src.diagnostic.plot import plotHeatMap
from src.parsers.bionano_file_reader import BionanoFileReader

matplotlib.use('TkAgg')

if __name__ == '__main__':
    for queryMeasurementError in ["100", "200", "500"]:
        queryName = f"query_meas{queryMeasurementError}"
        peaksFile = f"../.local_data/simulated/output/{queryName}_peaks.csv"
        simulationMetadataFilePath = f"../.local_data/simulated/{queryName}.sdata"

        simulationMetadata: DataFrame
        with open(simulationMetadataFilePath, "r") as simulationMetadataFile:
            simulationMetadata = BionanoFileReader(headersLinePrefix="#Fragment").readFile(
                simulationMetadataFile, ["ID", "Strand", "Start", "Stop"])

        peaks = pandas.read_csv(peaksFile)

        peaksWithSimulationMetadata = peaks.merge(simulationMetadata, left_on="queryId", right_on="ID")
        peaksWithSimulationMetadata["valid"] = peaksWithSimulationMetadata.apply(
            lambda row: ((row["Strand"] == "-") == row["reverseStrand"]) and int(row["Start"]) <= int(row["peakPosition"]) <= int(row["Stop"]), axis=1)

        validByResolutionAndBlur = peaksWithSimulationMetadata[["resolution", "blur", "valid"]].groupby(
            ["resolution", "blur"]).mean().unstack()

        x = validByResolutionAndBlur.columns.get_level_values("blur").unique()
        y = validByResolutionAndBlur.index.unique()
        plotHeatMap(
            validByResolutionAndBlur["valid"],
            f"../.local_data/simulated/output/{queryName}_heatmap.svg",
            x,
            y,
            title=f"Simulated molecules, {queryMeasurementError} measurement error")
