from tqdm import tqdm

from src.args import Args
from src.program import Program

if __name__ == '__main__':
    fileDirectory = "../.local_data/simulated/"
    referenceFile = fileDirectory + "chr1_r_unchanged.cmap"
    resolutionMultipliers = [0.5, 1, 2]
    blurs = [1, 2, 3]
    queryMeasurementErrors = [100, 200, 500, 1000]
    with tqdm(total=len(queryMeasurementErrors) * len(resolutionMultipliers) * len(blurs)) as progressBar:
        for measurementError in queryMeasurementErrors:
            queryName = f"query_meas{measurementError}"
            queryFile = fileDirectory + f"{queryName}.cmap"

            for resolutionMultiplier in resolutionMultipliers:
                resolution = int(measurementError * resolutionMultiplier)
                for blur in blurs:
                    outputFile = fileDirectory + f"output/secondary_{queryName}_r2_{resolution}_b2_{blur}.xmap"
                    args = Args.parse([
                        "-q", queryFile,
                        "-r", referenceFile,
                        "-o", outputFile,
                        "--primaryResolution", "1400",
                        "--primaryBlur", "1",
                        "--secondaryResolution", str(resolution),
                        "--secondaryBlur", str(blur),
                        "--peaksCount", "1"
                    ])

                    coma = Program(args)
                    coma.run()
                    progressBar.update()
