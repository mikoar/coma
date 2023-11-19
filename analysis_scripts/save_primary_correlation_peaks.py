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
    queryName = "chr1_simulated_1T_200T_scale001_meas100_res600_changed_8"
    referenceFile = "./.local_data/simulated/chr1_r_unchanged.cmap"
    queryFile = f"./.local_data/simulated/{queryName}.cmap"
    outputFile = f"./.local_data/simulated/output/{queryName}.xmap"
    peaksOutputFile = f"./.local_data/simulated/output/{queryName}_peaks.csv"

    with open(peaksOutputFile, "w") as f:
        f.write("resolution,blur,queryId,reverseStrand,score,peakPosition\n")

    resolutions = [230, 345, 460, 575, 690]
    blurs = [1, 2, 3]

    with tqdm(total=len(resolutions) * len(blurs), leave=False) as progressBar:
        for resolution in resolutions:
            for blur in blurs:
                args = Args.parse([
                    "-q", queryFile,
                    "-r", referenceFile,
                    "-o", outputFile,
                    "--primaryResolution", str(resolution),
                    "--primaryBlur", str(blur),
                    "--peaksCount", "1",
                    "--disableProgressBar"
                ])

                coma = Program(args, [PeaksCatcher(peaksOutputFile)])
                coma.run()
                progressBar.update()
