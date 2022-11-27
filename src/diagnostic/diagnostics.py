import os.path
from typing import TextIO

from src.diagnostic.alignment_plot import AlignmentPlot
from src.diagnostic.plot import plotCorrelation, plotRefinedCorrelation
from src.messaging.message_handler import MessageHandler
from src.messaging.messages import InitialAlignmentMessage, CorrelationResultMessage, AlignmentResultRowMessage


class DiagnosticsWriter:
    def __init__(self, outputFile: TextIO):
        self.outputDir = os.path.splitext(outputFile.name)[0] + "_diagnostics"
        os.makedirs(self.outputDir, exist_ok=True)

    def savePlot(self, fig, fileName: str):
        fig.savefig(os.path.join(self.outputDir, fileName), bbox_inches='tight',
                    pad_inches=0)


class PrimaryCorrelationDiagnosticsHandler(MessageHandler):
    messageType = InitialAlignmentMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: InitialAlignmentMessage):
        fig = plotCorrelation(message.data)
        self.writer.savePlot(fig, f"primary_cor_{message.data.query.moleculeId}.svg")


class SecondaryCorrelationDiagnosticsHandler(MessageHandler):
    messageType = CorrelationResultMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: CorrelationResultMessage):
        fig = plotRefinedCorrelation(message.initialAlignment, message.refinedAlignment)
        self.writer.savePlot(fig, f"secondary_cor{message.refinedAlignment.query.moleculeId}.svg")


class AlignmentPlotter(MessageHandler):
    messageType = AlignmentResultRowMessage

    def __init__(self, writer: DiagnosticsWriter):
        self.writer = writer

    def handle(self, message: AlignmentResultRowMessage):
        plot = AlignmentPlot(message.reference, message.query, message.alignment)
        self.writer.savePlot(plot.figure,
                             f"Alignment_ref_{message.reference.moleculeId}_query_{message.query.moleculeId}.svg")
