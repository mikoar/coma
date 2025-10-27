import os.path
import socket
from io import TextIOWrapper
from typing import List, TextIO, Iterable

import pandas as pd
from pandas import DataFrame, Series

from src.alignment.alignment_results import AlignmentResults
from src.args import Args
from src.correlation.bionano_alignment import BionanoAlignment
from src.parsers.bionano_file_reader import BionanoFileReader
from src.parsers.xmap_alignment_pair_parser import XmapAlignmentPairParser, BaseXmapAlignmentPairParser


class XmapReader:
    def __init__(self, pairParser: BaseXmapAlignmentPairParser = None) -> None:
        self.reader = BionanoFileReader()
        self.pairParser = pairParser or XmapAlignmentPairParser()

    def assign_xmap_ids(self, alignmentResults):
        for i, row in enumerate(alignmentResults.rows, start=1):
            setattr(row, "XmapEntryID", i)
            setattr(row, "xmapId", i)
            setattr(row, "alignmentId", i)

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

        return alignments.apply(self.rowParserFactory(), axis=1).tolist()

    def attach_alignment_strings(self, alignmentResults):
        for row in alignmentResults.rows:
            try:
                alignment_str = "".join(
                    f"({pair.reference.siteId},{pair.query.siteId})"
                    for pair in row.alignedPairs
                )
            except Exception:
                alignment_str = ""
            setattr(row, "Alignment", alignment_str)

    def writeAlignments(self, file, alignmentResults, args):

        self.assign_xmap_ids(alignmentResults)
        self.attach_alignment_strings(alignmentResults)

        file.write("# XMAP File Version: 0.2\n")
        file.write("#h XmapEntryID\tQryContigID\tRefContigID\tQryStartPos\tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment\n")
        file.write("#f int\tint\tint\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\tstring\tfloat\tfloat\tint\tstring\n")

        def _row_to_dict(row):
            q_start = float(row.queryEndPosition) if getattr(row, "orientation", "+") == "-" else float(row.queryStartPosition)
            q_end   = float(row.queryStartPosition) if getattr(row, "orientation", "+") == "-" else float(row.queryEndPosition)

            return {
                "XmapEntryID": int(getattr(row, "XmapEntryID")),
                "QryContigID": int(row.queryId),
                "RefContigID": int(row.referenceId),
                "QryStartPos": f"{q_start:.1f}",
                "QryEndPos": f"{q_end:.1f}",
                "RefStartPos": f"{float(row.referenceStartPosition):.1f}",
                "RefEndPos": f"{float(row.referenceEndPosition):.1f}",
                "Orientation": getattr(row, "orientation", "+"),
                "Confidence": f"{float(getattr(row, 'confidence', -1.0)):.2f}",
                "HitEnum": getattr(row, "cigarString", ""),
                "QryLen": f"{float(getattr(row, 'queryLength', 0.0)):.1f}",
                "RefLen": f"{float(getattr(row, 'referenceLength', 0.0)):.1f}",
                "LabelChannel": "1",
                "Alignment": getattr(row, "Alignment", ""),
            }

        for r in alignmentResults.rows:
            d = _row_to_dict(r)
            file.write(
                f"{d['XmapEntryID']}\t{d['QryContigID']}\t{d['RefContigID']}\t"
                f"{d['QryStartPos']}\t{d['QryEndPos']}\t{d['RefStartPos']}\t{d['RefEndPos']}\t"
                f"{d['Orientation']}\t{d['Confidence']}\t{d['HitEnum']}\t{d['QryLen']}\t{d['RefLen']}\t"
                f"{d['LabelChannel']}\t{d['Alignment']}\n"
            )

    def rowParserFactory(self):
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
    def argToString(arg):
        if isinstance(arg, TextIOWrapper):
            return arg.name
        if isinstance(arg, list):
            return " ".join(str(a) for a in arg)
        return str(arg)
