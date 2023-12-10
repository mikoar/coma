import os.path
from typing import TextIO, List

from matplotlib import pyplot as plt

from src.diagnostic.alignment_plot import AlignmentPlot
from src.diagnostic.plot import plotCorrelation, plotRefinedCorrelation
from src.extensions.extension import Extension
from src.extensions.messages import InitialAlignmentMessage, CorrelationResultMessage, AlignmentResultRowMessage, \
    MultipleAlignmentResultRowsMessage
from src.parsers.alignment_benchmark_reader import AlignmentBenchmarkReader


class DiagnosticsWriter:
    def __init__(self, outputFile: TextIO):
        if outputFile.name == '<stdout>':
            raise ValueError("outputFile is required")
        self.outputDir = os.path.splitext(outputFile.name)[0] + "_diagnostics"
        os.makedirs(self.outputDir, exist_ok=True)

    def savePlot(self, fig, fileName: str):
        fig.savefig(os.path.join(self.outputDir, fileName), bbox_inches='tight',
                    pad_inches=0)
        plt.close(fig)


class PrimaryCorrelationPlotter(Extension):
    messageType = InitialAlignmentMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: InitialAlignmentMessage):
        fig = plotCorrelation(message.data)
        self.writer.savePlot(fig, f"primary_cor_{message.data.query.moleculeId}"
                                  f"{'_reverse' if message.data.reverseStrand else ''}.svg")


class SecondaryCorrelationPlotter(Extension):
    messageType = CorrelationResultMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: CorrelationResultMessage):
        fig = plotRefinedCorrelation(message.initialAlignment, message.refinedAlignment)
        self.writer.savePlot(fig, f"secondary_cor_{message.refinedAlignment.query.moleculeId}_{message.index}.svg")


class AlignmentPlotter(Extension):
    messageType = AlignmentResultRowMessage

    def __init__(self, writer: DiagnosticsWriter, benchmarkReader: AlignmentBenchmarkReader, benchmarkFile: TextIO):
        self.writer = writer
        self.benchmarkReader = benchmarkReader
        self.benchmarkFile = benchmarkFile

    def handle(self, message: AlignmentResultRowMessage):
        if not message.alignment.alignedPairs:
            return

        benchmarkAlignment = self.getBenchmarkAlignment(message)
        plot = AlignmentPlot(message.reference, message.query, message.alignment, message.correlation,
                             benchmarkAlignment)

        self.writer.savePlot(plot.figure, f"Alignment_ref_{message.reference.moleculeId}_query"
                                          f"_{message.query.moleculeId}_{message.index}.svg")

    def getBenchmarkAlignment(self, message):
        if not self.benchmarkFile:
            return None
        return next(
            iter(self.benchmarkReader.read(self.benchmarkFile, queryIds=[message.query.moleculeId])))


class MultipleAlignmentsPlotter(Extension):
    messageType = MultipleAlignmentResultRowsMessage

    def __init__(self, writer: DiagnosticsWriter, benchmarkReader: AlignmentBenchmarkReader, benchmarkFile: TextIO):
        self.writer = writer
        self.benchmarkReader = benchmarkReader
        self.benchmarkAlignmentFile = benchmarkFile

    def handle(self, message: MultipleAlignmentResultRowsMessage):
        aligned = [m for m in message.messages if m.alignment.alignedPairs]
        if not aligned:
            return

        benchmarkAlignments = self.getBenchmarkAlignment([m.query.moleculeId for m in aligned])
        for m in aligned:
            plot = AlignmentPlot(m.reference, m.query, m.alignment, m.correlation,
                                 next((a for a in benchmarkAlignments if a.queryId == m.query.moleculeId), None))

            self.writer.savePlot(plot.figure, f"Alignment_ref_{m.reference.moleculeId}_query"
                                              f"_{m.query.moleculeId}_{m.index}.svg")

    def getBenchmarkAlignment(self, queryIds: List[int]):
        if not self.benchmarkAlignmentFile:
            return []
        return self.benchmarkReader.read(self.benchmarkAlignmentFile, queryIds=queryIds)
