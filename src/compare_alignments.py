import argparse
import sys
from typing import NamedTuple, TextIO, List

from src.diagnostic.alignment_comparer import AlignmentRowComparer, AlignmentComparer
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairWithDistanceParser
from src.parsers.xmap_reader import XmapReader


def main():
    parser = argparse.ArgumentParser(description="Compares optical map alignments.")
    parser.add_argument(dest="alignmentFiles", nargs=2, type=argparse.FileType("r"))
    parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"), required=True)
    parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"), required=True)
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"), default=sys.stdout)

    args: Args = parser.parse_args()  # type: ignore
    Program(args).run()


class Args(NamedTuple):
    alignmentFiles: List[TextIO]
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO


class Program:
    def __init__(self, args: Args):
        self.args = args
        self.sequenceReader = CmapReader()
        self.comparer = AlignmentComparer(AlignmentRowComparer())

    def run(self):
        references = self.sequenceReader.readReferences(self.args.referenceFile)
        queries = self.sequenceReader.readQueries(self.args.queryFile)
        alignmentReader = XmapReader(XmapAlignmentPairWithDistanceParser(references, queries))
        alignment1 = alignmentReader.readAlignments(self.args.alignmentFiles[0])
        alignment2 = alignmentReader.readAlignments(self.args.alignmentFiles[1])
        result = self.comparer.compare(alignment1, alignment2)
        result.write(self.args.outputFile)


if __name__ == '__main__':
    main()
