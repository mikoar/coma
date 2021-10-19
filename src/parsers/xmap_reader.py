from typing import List

import pandas as pd
from pandas import Series, DataFrame

from src.alignment.alignment_results import AlignmentResults
from src.correlation.alignment import Alignment
from src.parsers.bionano_file_reader import BionanoFileReader


class XmapReader:
    def __init__(self) -> None:
        self.reader = BionanoFileReader()

    def readAlignments(self, filePath: str, alignmentIds=None) -> List[Alignment]:
        alignments = self.reader.readFile(filePath,
                                          ["XmapEntryID", "QryContigID", "RefContigID", "QryStartPos",
                                           "QryEndPos", "RefStartPos", "RefEndPos", "Orientation",
                                           "Confidence", "QryLen", "Alignment"])
        if alignmentIds:
            alignments = alignments[alignments["XmapEntryID"].isin(alignmentIds)]

        return alignments.apply(self.__parseRow, axis=1).tolist()

    @staticmethod
    def writeAlignments(filePath: str, alignmentResults: AlignmentResults):
        with open(filePath, mode='w') as file:
            columns = {
                "#h": "#f",
                "XmapEntryID": "int",
                "QryContigID": "int",
                "RefContigID": "int",
                "RefStartPos": "float",
                "RefEndPos": "float",
                "Orientation": "string",
                "Confidence": "float",
                "Alignment": "string"
            }
            file.write("\t".join([columnName for columnName in columns.keys()]) + "\n")
            file.write("\t".join([columnType for columnType in columns.values()])+ "\n")

            dataFrame = DataFrame([{
                "QryContigID": row.queryId,
                "RefContigID": row.referenceId,
                "RefStartPos": "{:.1f}".format(row.referenceStartPosition),
                "RefEndPos": "{:.1f}".format(row.referenceEndPosition),
                "Orientation": "-" if row.reverseStrand else "+",
                "Confidence": "{:.2f}".format(row.score),
                "Alignment": "".join(
                    [f"({pair.referencePositionIndex},{pair.queryPositionIndex})" for pair in row.alignedPairs]),
            } for row in alignmentResults.rows], index=pd.RangeIndex(start=1, stop=len(alignmentResults.rows) + 1))

            dataFrame.to_csv(file, sep='\t', header=False, mode="a", line_terminator="\n")

    @staticmethod
    def __parseRow(row: Series):
        return Alignment.parse(row["XmapEntryID"], row["QryContigID"], row["RefContigID"], row["QryStartPos"],
                               row["QryEndPos"], row["RefStartPos"], row["RefEndPos"], row["Orientation"],
                               row["Confidence"], row["QryLen"], row["Alignment"])
