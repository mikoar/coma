from __future__ import annotations

import argparse
import os
from typing import NamedTuple, TextIO, List

from tqdm import tqdm

from src.diagnostic.alignment_plot import BenchmarkAlignmentPlot
from src.diagnostic.diagnostics import DiagnosticsWriter
from src.parsers.alignment_benchmark_reader import AlignmentBenchmarkReader
from src.parsers.cmap_reader import CmapReader
from src.parsers.simulation_alignment_pair_parser import SimulationAlignmentPairWithDistanceParser
from src.parsers.simulation_data_as_xmap_reader import SimulationDataAsXmapReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairWithDistanceParser
from src.parsers.xmap_reader import XmapReader


def main():
    parser = argparse.ArgumentParser(description="Plots optical map alignments.")
    parser.add_argument(dest="alignmentFile", nargs=1, type=argparse.FileType("r"))
    parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"), required=True)
    parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"), required=True)
    parser.add_argument("-o", "--output", dest="outputFile", type=argparse.FileType("w"), required=True)
    parser.add_argument("-qId", "--queryIDs", dest="queryIds", type=int, nargs="*")
    parser.add_argument("-n", "--maxCount", dest="maxCount", type=int, default=None)

    args: Args = parser.parse_args()  # type: ignore
    Program(args).run()


class Args(NamedTuple):
    alignmentFile: List[TextIO]
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO
    queryIds: List[int]
    maxCount: int | None


class Program:
    def __init__(self, args: Args):
        self.args = args
        self.sequenceReader = CmapReader()
        self.writer = DiagnosticsWriter(args.outputFile)

    def run(self):
        references = self.sequenceReader.readReferences(self.args.referenceFile)
        queries = self.sequenceReader.readQueries(self.args.queryFile)
        pairParser = XmapAlignmentPairWithDistanceParser(references, queries)
        benchmarkReader = AlignmentBenchmarkReader(
            XmapReader(pairParser),
            SimulationDataAsXmapReader(SimulationAlignmentPairWithDistanceParser(references, queries)))
        alignments = benchmarkReader.read(self.args.alignmentFile[0], self.args.queryIds)
        for alignment in tqdm(alignments[:self.args.maxCount or len(alignments)]):
            reference = next(r for r in references if r.moleculeId == alignment.referenceId)
            query = next(q for q in queries if q.moleculeId == alignment.queryId)
            plot = BenchmarkAlignmentPlot(
                reference,
                query,
                alignment)
            self.writer.savePlot(plot.figure, f"r{reference.moleculeId}_q{query.moleculeId}_{os.path.basename(self.args.outputFile.name)}")


if __name__ == '__main__':
    main()
