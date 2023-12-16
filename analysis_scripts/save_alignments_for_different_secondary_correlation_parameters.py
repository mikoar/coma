from tqdm import tqdm

from src.args import Args as ComaArgs
from src.compare_alignments import Program as CompareAlignments, Args as CompareArgs
from src.program import Program as Coma

if __name__ == '__main__':
    fileDirectory = "../.local_data/simulated/"
    referenceFile = fileDirectory + "chr1_r_unchanged.cmap"
    resolutions = [50, 100, 200, 400, 800]
    blurs = [1, 2, 3, 4, 5]
    with tqdm(total=len(resolutions) * len(blurs)) as progressBar:
        # for measurementError in queryMeasurementErrors:
        queryName = f"filtered_molecules_1500_200T_scale001_meas500_res600_fpr000002"
        queryFile = fileDirectory + f"{queryName}.cmap"

        for resolution in resolutions:
            for blur in blurs:
                outputFileName = fileDirectory + f"output/2cc_params/{queryName}_r2_{resolution}_b2_{blur}"
                xmapFile = outputFileName + ".xmap"
                comaArgs = ComaArgs.parse([
                    "-q", queryFile,
                    "-r", referenceFile,
                    "-o", xmapFile,
                    "--primaryResolution", "1400",
                    "--primaryBlur", "1",
                    "--secondaryResolution", str(resolution),
                    "--secondaryBlur", str(blur),
                    "--peaksCount", "3",
                    "--maxPairDistance", "1000",
                    "--distancePenaltyMultiplier", "0.5"
                ])

                coma = Coma(comaArgs)
                coma.run()

                compareArgs = CompareArgs.parse([
                    xmapFile,
                    fileDirectory + f"{queryName}.sdata",
                    "-o", outputFileName + ".tsv",
                    "-q", queryFile,
                    "-r", referenceFile,
                    "-d"
                ])
                compare = CompareAlignments(compareArgs)
                compare.run()
                progressBar.update()
