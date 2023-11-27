from typing import List, TextIO, Iterable

from pandas import Series

from src.correlation.simulated_alignment import SimulatedAlignment
from src.diagnostic.xmap_alignment import XmapAlignment, XmapAlignedPair
from src.parsers.bionano_file_reader import BionanoFileReader


class SimulationDataAsXmapReader:
    def __init__(self, reader: BionanoFileReader = None) -> None:
        self.reader = reader or BionanoFileReader(headersLinePrefix="#Fragment")

    def readAlignments(self, file: TextIO, queryIds: Iterable[int] = None) -> List[XmapAlignment]:
        alignments = self.reader.readFile(
            file,
            ["ID", "Reference", "Strand", "Start", "Stop", "SimuInfoDetail", "Size"])

        if queryIds:
            alignments = alignments[alignments["ID"].isin(queryIds)]

        return alignments.apply(self.__rowParserFactory(), axis=1).tolist()

    def __rowParserFactory(self):
        def parseRow(row: Series):
            queryId = int(row["ID"])
            referenceId = int(row["Reference"])
            reverseStrand = row["Strand"] == "-"
            return SimulatedAlignment.parse(queryId, queryId, referenceId, 0,
                                            row["Size"], row["Start"], row["Stop"], reverseStrand,
                                            9999., "", row["Size"], row["Size"],
                                            self.__parsePairs(row["SimuInfoDetail"], reverseStrand))

        return parseRow

    @staticmethod
    def __parsePairs(simulationDetail: str, reverseStrand: bool) -> List[XmapAlignedPair]:
        pairs = [pair for pairs in
                 [SimulationDataAsXmapReader.__parsePair(x, i + 1) for i, x in enumerate(simulationDetail.split(";")) if x != "FP"]
                 for pair in pairs]
        if reverseStrand:
            pairs.reverse()
        return pairs

    @staticmethod
    def __parsePair(positionString: str, queryId: int):
        if "," in positionString:
            splitPositionStrings = positionString.split(",")
            return [XmapAlignedPair.create(s.split(":")[1], str(queryId)) for s in splitPositionStrings if s != "FP"]
        else:
            return [XmapAlignedPair.create(positionString.split(":")[1], str(queryId))]
