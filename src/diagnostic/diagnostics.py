import os.path
from abc import ABC, abstractmethod
from typing import TextIO

from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.optical_map import InitialAlignment, CorrelationResult
from src.correlation.plot import plotCorrelation


class Diagnostics(ABC):
    @staticmethod
    def create(diagnosticsEnabled: bool, benchmarkAlignmentFile: TextIO, outputFile: TextIO):
        return _DiagnosticsImpl(benchmarkAlignmentFile, outputFile) if diagnosticsEnabled else _NullDiagnostics()

    @abstractmethod
    def primaryCorrelation(self, primaryCorrelation: InitialAlignment):
        pass

    @abstractmethod
    def secondaryCorrelation(self, secondaryCorrelation: CorrelationResult):
        pass

    @abstractmethod
    def alignment(self, alignment: AlignmentResultRow):
        pass


class _DiagnosticsImpl(Diagnostics):
    def __init__(self, benchmarkAlignmentFile: TextIO, outputFile: TextIO):
        self.benchmarkAlignmentFile = benchmarkAlignmentFile
        self.outputDir = os.path.join(os.path.dirname(outputFile.name), "diagnostics")
        os.makedirs(self.outputDir, exist_ok=True)

    def primaryCorrelation(self, primaryCorrelation: InitialAlignment):
        fig = plotCorrelation(primaryCorrelation)
        self.__savePlot(fig, f"primary_cor_{primaryCorrelation.query.moleculeId}.svg")

    def secondaryCorrelation(self, secondaryCorrelation: CorrelationResult):
        fig = plotCorrelation(secondaryCorrelation)
        self.__savePlot(fig, f"secondary_cor_{secondaryCorrelation.query.moleculeId}.svg")

    def alignment(self, alignment: AlignmentResultRow):
        pass

    def __savePlot(self, fig, fileName: str):
        fig.savefig(os.path.join(self.outputDir, fileName), bbox_inches='tight',
                    pad_inches=0)


class _NullDiagnostics(Diagnostics):
    def secondaryCorrelation(self, secondaryCorrelation: CorrelationResult):
        pass

    def alignment(self, alignment: AlignmentResultRow):
        pass

    def primaryCorrelation(self, primaryCorrelation: InitialAlignment):
        pass
