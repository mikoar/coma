from typing import List, TextIO, Iterable

from pandas import Series

from src.correlation.simulated_alignment import SimulatedAlignment
from src.diagnostic.benchmark_alignment import BenchmarkAlignment
from src.parsers.bionano_file_reader import BionanoFileReader
from src.parsers.simulation_alignment_pair_parser import BaseSimulationAlignmentPairParser, \
    SimulationAlignmentPairParser


class SimulationDataAsXmapReader:
    def __init__(self, pairParser: BaseSimulationAlignmentPairParser = None, reader: BionanoFileReader = None) -> None:
        self.pairParser = pairParser or SimulationAlignmentPairParser()
        self.reader = reader or BionanoFileReader(headersLinePrefix="#Fragment")

    def readAlignments(self, file: TextIO, queryIds: Iterable[int] = None) -> List[BenchmarkAlignment]:
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
                                            self.pairParser.parse(row["SimuInfoDetail"], queryId, referenceId, reverseStrand))

        return parseRow
