import os.path
from typing import TextIO

from matplotlib import pyplot as plt

from src.diagnostic.alignment_plot import AlignmentPlot
from src.diagnostic.plot import plotCorrelation, plotRefinedCorrelation
from src.messaging.message_handler import MessageHandler
from src.messaging.messages import InitialAlignmentMessage, CorrelationResultMessage, AlignmentResultRowMessage
from src.parsers.xmap_reader import XmapReader


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


class PrimaryCorrelationDiagnosticsHandler(MessageHandler):
    messageType = InitialAlignmentMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: InitialAlignmentMessage):
        fig = plotCorrelation(message.data)
        self.writer.savePlot(fig, f"primary_cor_{message.data.query.moleculeId}"
                                  f"{'_reverse' if message.data.reverseStrand else ''}.svg")


class SecondaryCorrelationDiagnosticsHandler(MessageHandler):
    messageType = CorrelationResultMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: CorrelationResultMessage):
        fig = plotRefinedCorrelation(message.initialAlignment, message.refinedAlignment)
        self.writer.savePlot(fig, f"secondary_cor{message.refinedAlignment.query.moleculeId}_{message.index}.svg")


class AlignmentPlotter(MessageHandler):
    messageType = AlignmentResultRowMessage

    def __init__(self, writer: DiagnosticsWriter, xmapReader: XmapReader, benchmarkAlignmentFile: TextIO):
        self.writer = writer
        self.xmapReader = xmapReader
        self.benchmarkAlignmentFile = benchmarkAlignmentFile

    def handle(self, message: AlignmentResultRowMessage):
        if not message.alignment.alignedPairs:
            return

        benchmarkAlignment = self.getBenchmarkAlignment(message)
        plot = AlignmentPlot(message.reference, message.query, message.alignment, message.correlation,
                             benchmarkAlignment)

        self.writer.savePlot(plot.figure, f"Alignment_ref_{message.reference.moleculeId}_query"
                                          f"_{message.query.moleculeId}_{message.index}.svg")

    def getBenchmarkAlignment(self, message):
        if not self.benchmarkAlignmentFile:
            return None
        return next(
            iter(self.xmapReader.readAlignments(self.benchmarkAlignmentFile, queryIds=[message.query.moleculeId])))
