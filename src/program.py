from __future__ import annotations

import sys
from typing import List

from src.alignment.alignment_results import AlignmentResults
from src.application_service import ApplicationServiceFactory
from src.args import Args
from src.diagnostic.diagnostics import DiagnosticsWriter, PrimaryCorrelationDiagnosticsHandler, \
    SecondaryCorrelationDiagnosticsHandler, AlignmentPlotter, MultipleAlignmentsPlotter
from src.messaging.dispatcher import Dispatcher
from src.messaging.message_handler import MessageHandler
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairWithDistanceParser
from src.parsers.xmap_reader import XmapReader


def main():
    args = Args.parse()
    Program(args).run()


class Program:
    def __init__(self, args: Args, extensions: List[MessageHandler] = None):
        self.args = args
        self.__readMaps()
        self.xmapReader = XmapReader(XmapAlignmentPairWithDistanceParser(self.referenceMaps, self.queryMaps))
        self.dispatcher = Dispatcher(extensions)
        self.applicationService = ApplicationServiceFactory().create(args, self.dispatcher)
        if args.diagnosticsEnabled:
            writer = DiagnosticsWriter(args.outputFile)
            self.dispatcher.addHandler(PrimaryCorrelationDiagnosticsHandler(writer))
            self.dispatcher.addHandler(SecondaryCorrelationDiagnosticsHandler(writer))
            self.dispatcher.addHandler(AlignmentPlotter(writer, self.xmapReader, args.benchmarkAlignmentFile))
            self.dispatcher.addHandler(MultipleAlignmentsPlotter(writer, self.xmapReader, args.benchmarkAlignmentFile))

    def run(self):
        alignmentResultRows = self.applicationService.execute(self.referenceMaps, self.queryMaps)
        alignmentResult = AlignmentResults.create(self.args.referenceFile.name, self.args.queryFile.name,
                                                  alignmentResultRows)
        self.xmapReader.writeAlignments(self.args.outputFile, alignmentResult)
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
