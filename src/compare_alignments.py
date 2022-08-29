import argparse
import sys
from typing import NamedTuple, TextIO

from src.alignment.alignment_comparer import AlignmentRowComparer, AlignmentComparer
from src.parsers.xmap_reader import XmapReader


def main():
    parser = argparse.ArgumentParser(description="Compares optical map alignments.")
    parser.add_argument("-r", "--reference", dest="referenceFile", type=argparse.FileType("r"))
    parser.add_argument("-q", "--query", dest="queryFile", type=argparse.FileType("r"))
    parser.add_argument("-o", "--output", dest="outputFile", nargs="?", type=argparse.FileType("w"), default=sys.stdout)

    args: Args = parser.parse_args()  # type: ignore
    Program(args).run()


class Args(NamedTuple):
    referenceFile: TextIO
    queryFile: TextIO
    outputFile: TextIO


class Program:
    def __init__(self, args: Args):
        self.args = args
        self.reader = XmapReader()
        self.comparer = AlignmentComparer(AlignmentRowComparer())

    def run(self):
        reference = self.reader.readAlignments(self.args.referenceFile)
        query = self.reader.readAlignments(self.args.queryFile)
        result = self.comparer.compare(reference, query)
        result.write(self.args.outputFile)


if __name__ == '__main__':
    main()
