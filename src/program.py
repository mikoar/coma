from __future__ import annotations

import sys
from typing import List

from src.alignment.alignment_results import AlignmentResults
from src.args import Args
from src.diagnostic.diagnostics import DiagnosticsWriter, PrimaryCorrelationPlotter, \
    SecondaryCorrelationPlotter, AlignmentPlotter, MultipleAlignmentsPlotter
from src.extensions.dispatcher import Dispatcher
from src.extensions.extension import Extension
from src.parsers.alignment_benchmark_reader import AlignmentBenchmarkReader
from src.parsers.cmap_reader import CmapReader
from src.parsers.simulation_alignment_pair_parser import SimulationAlignmentPairWithDistanceParser
from src.parsers.simulation_data_as_xmap_reader import SimulationDataAsXmapReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairWithDistanceParser
from src.parsers.xmap_reader import XmapReader
from src.workflow_coordinator_factory import WorkflowCoordinatorFactory


def main():
    args = Args.parse()
    Program(args).run()


class Program:
    def __init__(self, args: Args, extensions: List[Extension] = None):
        self.args = args
        self.__readMaps()
        self.xmapReader = XmapReader(XmapAlignmentPairWithDistanceParser(self.referenceMaps, self.queryMaps))
        self.dispatcher = Dispatcher(extensions)
        self.workflowCoordinator = WorkflowCoordinatorFactory(args, self.dispatcher, self.xmapReader).create()
        if args.diagnosticsEnabled:
            writer = DiagnosticsWriter(args.outputFile)
            simulationDataReader = SimulationDataAsXmapReader(
                SimulationAlignmentPairWithDistanceParser(self.referenceMaps, self.queryMaps))
            benchmarkReader = AlignmentBenchmarkReader(self.xmapReader, simulationDataReader)
            self.dispatcher.addExtension(PrimaryCorrelationPlotter(writer))
            self.dispatcher.addExtension(SecondaryCorrelationPlotter(writer))
            self.dispatcher.addExtension(AlignmentPlotter(writer, benchmarkReader, args.benchmarkAlignmentFile))
            self.dispatcher.addExtension(
                MultipleAlignmentsPlotter(writer, benchmarkReader, args.benchmarkAlignmentFile))

    def run(self):
        alignmentResultRows = self.workflowCoordinator.execute(self.referenceMaps, self.queryMaps)
        alignmentResult = AlignmentResults.create(self.args.referenceFile.name, self.args.queryFile.name,
                                                  alignmentResultRows)
        self.xmapReader.writeAlignments(self.args.outputFile, alignmentResult, self.args)
        if self.args.outputFile is not sys.stdout:
            self.args.outputFile.close()
        return alignmentResult

    def __readMaps(self):
        cmapReader = CmapReader()
        with self.args.referenceFile:
            self.referenceMaps = cmapReader.readReferences(self.args.referenceFile, self.args.referenceIds)
        with self.args.queryFile:
            self.queryMaps = list(
                map(lambda q: q.trim(), cmapReader.readQueries(self.args.queryFile, self.args.queryIds)))


if __name__ == '__main__':
    main()
