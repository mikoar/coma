from tqdm import tqdm

from src.args import Args
from src.extensions.extension import Extension
from src.extensions.messages import CorrelationResultMessage
from src.program import Program


class PeaksCatcher(Extension):
    messageType = CorrelationResultMessage

    def __init__(self, filePath: str):
        self.filePath = filePath

    def handle(self, message: CorrelationResultMessage):
        peak = message.initialAlignment.maxPeak
        if peak:
            with open(self.filePath, "a") as f:
                f.write(f"{message.initialAlignment.resolution},"
                        f"{message.initialAlignment.blur},"
                        f"{message.initialAlignment.query.moleculeId},"
                        f"{message.initialAlignment.reverseStrand},"
                        f"{peak.score:.2f},"
                        f"{peak.position}\n")


if __name__ == '__main__':
    fileDirectory = "../.local_data/simulated/"
    referenceFile = fileDirectory + "chr1_r_unchanged.cmap"
    resolutions = [500, 700, 1000, 1400, 2000]
    blurs = [1, 2, 3, 4]
    queryMeasurementErrors = ["100", "200", "500"]
    with tqdm(total=len(queryMeasurementErrors) * len(resolutions) * len(blurs)) as progressBar:
        for measurementError in queryMeasurementErrors:
            queryName = f"query_meas{measurementError}"
            queryFile = fileDirectory + f"{queryName}.cmap"
            peaksOutputFile = fileDirectory + f"{queryName}_peaks.csv"

            with open(peaksOutputFile, "w") as f:
                f.write("resolution,blur,queryId,reverseStrand,score,peakPosition\n")
                for resolution in resolutions:
                    for blur in blurs:
                        outputFile = fileDirectory + f"{queryName}_r1_{resolution}_b1{blur}.xmap"
                        args = Args.parse([
                            "-q", queryFile,
                            "-r", referenceFile,
                            "-o", outputFile,
                            "--primaryResolution", str(resolution),
                            "--primaryBlur", str(blur),
                            "--peaksCount", "1"
                        ])

                        coma = Program(args, [PeaksCatcher(peaksOutputFile)])
                        coma.run()
                        progressBar.update()
