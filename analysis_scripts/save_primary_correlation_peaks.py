from tqdm import tqdm


from src.args import Args
from src.extensions.extension import Extension
from src.extensions.messages import CorrelationResultMessage, AlignmentResultRowMessage, MultipleAlignmentResultRowsMessage
from src.program import Program


class PeaksCatcher(Extension):
    messageType = CorrelationResultMessage

    def __init__(self, filePath: str):
        self.filePath = filePath

    def handle(self, message: CorrelationResultMessage):
        peak = message.initialAlignment.maxPeak
        if peak:
            with open(self.filePath, "a") as f:
                f.write(f"{message.initialAlignment.resolution};"
                        f"{message.initialAlignment.blur};"
                        f"{message.initialAlignment.query.moleculeId};"
                        f"{message.initialAlignment.reverseStrand};"
                        f"{peak.score:.2f};"
                        f"{peak.position}\n")

class SegmentsCatcher(Extension):
    messageType = AlignmentResultRowMessage

    def __init__(self, filePath: str):
        self.filePath = filePath

    def handle(self, message: AlignmentResultRowMessage):
        peak = message.correlation.maxPeak
        if peak:
            with open(self.filePath, "a") as f:
                f.write(
                        f"{message.correlation.query.moleculeId};"
                        f"{message.correlation.reverseStrand};"
                        f"{len(message.alignment.segments)};"
                        f"{message.alignment.segments};"
                        f"{message.alignment.queryStartPosition};"
                        f"{message.alignment.queryEndPosition};"
                        f"{message.alignment.referenceStartPosition};"
                        f"{message.alignment.referenceEndPosition};"
                        f"{peak.score:.2f};"
                        f"{peak.position}\n")
                
class MultipleSegmentsCatcher(Extension):
    messageType = MultipleAlignmentResultRowsMessage

    def __init__(self, filePath: str):
        self.filePath = filePath

    def handle(self, message: AlignmentResultRowMessage):
        peak = message.correlation.maxPeak
        if peak:
            with open(self.filePath, "a") as f:
                f.write(
                        f"{message.correlation.query.moleculeId};"
                        f"{message.correlation.reverseStrand};"
                        f"{len(message.alignment.segments)};"
                        f"{message.alignment.segments};"
                        f"{message.alignment.queryStartPosition};"
                        f"{message.alignment.queryEndPosition};"
                        f"{message.alignment.referenceStartPosition};"
                        f"{message.alignment.referenceEndPosition};"
                        f"{peak.score:.2f};"
                        f"{peak.position}\n")

if __name__ == '__main__':
    fileDirectory = "../OM_data/"
    referenceFile = fileDirectory + "chr1_r_new_fandom.cmap"
    thresholds = [5, 10, 30]
    peaksNumbers = [1, 3, 5]
    with tqdm(total=len(thresholds) * len(peaksNumbers)) as progressBar:
        queryFile = fileDirectory + "query.cmap"
        for threshold in thresholds:
            for peak in peaksNumbers:
                segmentsOutputFile = fileDirectory + f"output/query_peaks_th_{threshold}_pc_{peak}.csv"
                with open(segmentsOutputFile, "w") as f:
                    f.write("queryId;reverseStrand;segmentNb;segment;segmentQueryStart;segmentQueryEnd;segmentReferenceStart;segmentReferenceEnd;score;peakPosition\n")

                outputFile = fileDirectory + f"output/query_peaks_th_{threshold}_pc_{peak}.xmap"
                args = Args.parse([
                    "-q", queryFile,
                    "-r", referenceFile,
                    "-o", outputFile,
                    "--peakHeightThreshold", str(threshold),
                    "--peaksCount", str(peak)
                ])
                coma = Program(args, [SegmentsCatcher(segmentsOutputFile)])
                coma.run()
                progressBar.update()
