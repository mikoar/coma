import os.path
from typing import TextIO

from src.diagnostic.plot import plotCorrelation, plotRefinedCorrelation
from src.messaging.message_handler import MessageHandler
from src.messaging.messages import InitialAlignmentMessage, CorrelationResultMessage


class DiagnosticsPlotter:
    def __init__(self, outputFile: TextIO):
        self.outputDir = os.path.splitext(outputFile.name)[0] + "_diagnostics"
        os.makedirs(self.outputDir, exist_ok=True)

    def savePlot(self, fig, fileName: str):
        fig.savefig(os.path.join(self.outputDir, fileName), bbox_inches='tight',
                    pad_inches=0)


class PrimaryCorrelationDiagnosticsHandler(MessageHandler):
    messageType = InitialAlignmentMessage

    def __init__(self, plotter: DiagnosticsPlotter):
        self.writer = plotter

    def handle(self, message: InitialAlignmentMessage):
        fig = plotCorrelation(message.data)
        self.writer.savePlot(fig, f"primary_cor_{message.data.query.moleculeId}.svg")


class SecondaryCorrelationDiagnosticsHandler(MessageHandler):
    messageType = CorrelationResultMessage

    def __init__(self, plotter: DiagnosticsPlotter):
        self.writer = plotter

    def handle(self, message: CorrelationResultMessage):
        fig = plotRefinedCorrelation(message.initialAlignment, message.refinedAlignment)
        self.writer.savePlot(fig, f"secondary_cor{message.refinedAlignment.query.moleculeId}.svg")
