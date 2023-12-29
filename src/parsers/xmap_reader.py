import os.path
import socket
from io import TextIOWrapper
from typing import List, TextIO, Iterable

import pandas as pd
from pandas import DataFrame, Series

from src.args import Args
from src.alignment.alignment_results import AlignmentResults
from src.correlation.bionano_alignment import BionanoAlignment
from src.parsers.bionano_file_reader import BionanoFileReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairParser, BaseXmapAlignmentPairParser


class XmapReader:
    def __init__(self, pairParser: BaseXmapAlignmentPairParser = None) -> None:
        self.reader = BionanoFileReader()
        self.pairParser = pairParser or XmapAlignmentPairParser()

    def readAlignments(self, file: TextIO, alignmentIds: Iterable[int] = None, queryIds: Iterable[int] = None) -> \
            List[BionanoAlignment]:
        alignments = self.reader.readFile(file,
                                          ["XmapEntryID", "QryContigID", "RefContigID", "QryStartPos",
                                           "QryEndPos", "RefStartPos", "RefEndPos", "Orientation",
                                           "Confidence", "HitEnum", "QryLen", "RefLen", "Alignment"])
        if alignmentIds:
            alignments = alignments[alignments["XmapEntryID"].isin(alignmentIds)]

        if queryIds:
            alignments = alignments[alignments["QryContigID"].isin(queryIds)]

        return alignments.apply(self.__rowParserFactory(), axis=1).tolist()

    def writeAlignments(self, file: TextIO, alignmentResults: AlignmentResults, args: Args):
        columns = {
            "#h": "#f",
            "XmapEntryID": "int",
            "QryContigID": "int",
            "RefContigID": "int",
            "QryStartPos": "float",
            "QryEndPos": "float",
            "RefStartPos": "float",
            "RefEndPos": "float",
            "Orientation": "string",
            "Confidence": "float",
            "HitEnum": "string",
            "QryLen": "float",
            "RefLen": "float",
            "AlignedRest" : "string",
            "LabelChannel": "int",
            "Alignment": "string"
        }
        file.writelines(line + "\n" for line in [
            f"# hostname={socket.gethostname()}",
            "# coma " + " ".join([f"--{k} {self.__argToString(v)}" for k, v in vars(args).items()]),
            "# XMAP File Version:\t0.2",
            f"# Reference Maps From:\t{os.path.abspath(alignmentResults.referenceFilePath)}",
            f"# Query Maps From:\t{os.path.abspath(alignmentResults.queryFilePath)}",
            "\t".join([columnName for columnName in columns.keys()]),
            "\t".join([columnType for columnType in columns.values()])])

        dataFrame = DataFrame([{
            "QryContigID": row.queryId,
            "RefContigID": row.referenceId,
            "QryStartPos": "{:.1f}".format(row.queryStartPosition),
            "QryEndPos": "{:.1f}".format(row.queryEndPosition),
            "RefStartPos": "{:.1f}".format(row.referenceStartPosition),
            "RefEndPos": "{:.1f}".format(row.referenceEndPosition),
            "Orientation": row.orientation,
            "Confidence": "{:.2f}".format(row.confidence),
            "HitEnum": row.cigarString,
            "QryLen": "{:.1f}".format(row.queryLength),
            "RefLen": "{:.1f}".format(row.referenceLength),
            "AlignedRest" : "{}".format(row.alignedRest),
            "LabelChannel": 1,
            "Alignment": "".join(
                [f"({pair.reference.siteId},{pair.query.siteId})" for pair in row.alignedPairs]),
        } for row in alignmentResults.rows], index=pd.RangeIndex(start=1, stop=len(alignmentResults.rows) + 1))
        dataFrame.to_csv(file, sep='\t', header=False, mode="a", line_terminator='\n')

    def __rowParserFactory(self):
        def parseRow(row: Series):
            queryId = int(row["QryContigID"])
            referenceId = int(row["RefContigID"])
            reverseStrand = row["Orientation"] == "-"
            return BionanoAlignment.parse(row["XmapEntryID"], queryId, referenceId, row["QryStartPos"],
                                          row["QryEndPos"], row["RefStartPos"], row["RefEndPos"], reverseStrand,
                                          row["Confidence"], row["HitEnum"], row["QryLen"], row["RefLen"],
                                          self.pairParser.parse(row["Alignment"], queryId, referenceId, reverseStrand))

        return parseRow

    @staticmethod
    def __argToString(arg):
        if isinstance(arg, TextIOWrapper):
            return arg.name
        if isinstance(arg, list):
            return " ".join(str(a) for a in arg)
        return str(arg)
