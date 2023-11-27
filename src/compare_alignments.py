import argparse
import itertools
import sys
from typing import NamedTuple, TextIO, List

from src.diagnostic.alignment_comparer import AlignmentRowComparer, AlignmentComparer
from src.parsers.cmap_reader import CmapReader
from src.parsers.simulation_data_as_xmap_reader import SimulationDataAsXmapReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairWithDistanceParser, XmapAlignmentPairParser
from src.parsers.xmap_reader import XmapReader


def main():
    parser = argparse.ArgumentParser(description="Compares optical map alignments.")
    parser.add_argument(dest="alignmentFiles", nargs=2, type=argparse.FileType("r"))
    parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"), required=True)
    parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"), required=True)
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"), default=sys.stdout)
    parser.add_argument("-d", "--includePositions", dest="includePositions", action="store_true",
                        help="Additionally writes genomic positions for each label in the output.")

    args: Args = parser.parse_args()  # type: ignore
    Program(args).run()


class Args(NamedTuple):
    alignmentFiles: List[TextIO]
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO
    includePositions: bool


class Program:
    def __init__(self, args: Args):
        self.args = args
        self.sequenceReader = CmapReader()
        self.comparer = AlignmentComparer(AlignmentRowComparer())

    def run(self):
        if self.args.includePositions:
            references = self.sequenceReader.readReferences(self.args.referenceFile)
            queries = self.sequenceReader.readQueries(self.args.queryFile)
            pairParser = XmapAlignmentPairWithDistanceParser(references, queries)
        else:
            pairParser = XmapAlignmentPairParser()
        alignmentReader = XmapReader(pairParser)
        simulatedAlignmentReader = SimulationDataAsXmapReader()
        alignment1 = self.__read(self.args.alignmentFiles[0], alignmentReader, simulatedAlignmentReader)
        alignment2 = self.__read(self.args.alignmentFiles[1], alignmentReader, simulatedAlignmentReader)
        result = self.comparer.compare(alignment1, alignment2)
        result.write(self.args.outputFile)

    @staticmethod
    def __read(file: TextIO, alignmentReader: XmapReader, simulatedAlignmentReader: SimulationDataAsXmapReader):
        headers = "\n".join(itertools.islice(itertools.takewhile(lambda line: line.startswith("#"), file), 10))
        file.seek(0, 0)
        if "XMAP" in headers:
            return alignmentReader.readAlignments(file)
        if "SimuInfoDetail" in headers:
            return simulatedAlignmentReader.readAlignments(file)
        else:
            raise Exception(f"File {file.name} is in unknown format. Either XMAP or SDATA formats are supported.")


if __name__ == '__main__':
    main()
