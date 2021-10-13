import itertools
import re
from typing import List

import pandas
from pandas.core.frame import DataFrame
from pandas.core.series import Series

from src.correlation.alignment import Alignment
from src.correlation.optical_map import VectorisedOpticalMap
from src.correlation.sequence_generator import SequenceGenerator


class BionanoFileReader:
    def readFile(self, filePath, columns: List[str]) -> DataFrame:
        return pandas.read_csv(
            filePath,
            comment="#",
            delimiter="\t",
            names=self.__getColumnNames(filePath),
            usecols=columns)  # type: ignore

    def __getColumnNames(self, filePath):
        with open(filePath) as file:
            gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
            header_line = list(itertools.islice(gen, 1))[0].strip()
            names = re.split(r'\s+', header_line)[1:]
        return names


class CmapReader:
    def __init__(self, sequenceGenerator: SequenceGenerator) -> None:
        self.sequenceGenerator = sequenceGenerator
        self.reader = BionanoFileReader()

    def readQueries(self, filePath: str, moleculeIds: List[int] = []):
        return self.__read(filePath, moleculeIds)

    def readQuery(self, filePath: str, moleculeId: int):
        return self.__read(filePath, [moleculeId])[0]

    def readReference(self, filePath: str, chromosome: int = 1):
        return self.__read(filePath, [chromosome])[0]

    def readReferences(self, filePath: str, chromosomes: List[int] = []):
        return self.__read(filePath, chromosomes)

    def __read(self, filePath, moleculeIds=None) -> List[VectorisedOpticalMap]:
        maps = self.reader.readFile(filePath, ["CMapId", "Position", "ContigLength"])

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup, self.sequenceGenerator)
        return opticalMaps.tolist()

    @staticmethod
    def __parseCmapRowsGroup(group, sequenceGenerator: SequenceGenerator):
        moleculeId = group["CMapId"].iloc[0]
        length = int(group["ContigLength"].iloc[0])
        positions = group["Position"].sort_values().tolist()
        sequence = sequenceGenerator.positionsToSequence(positions)
        return VectorisedOpticalMap(moleculeId, length, positions, sequence, sequenceGenerator.resolution)


class AlignmentReader:
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
    def __parseRow(row: Series):
        return Alignment.parse(row["XmapEntryID"], row["QryContigID"], row["RefContigID"], row["QryStartPos"],
                               row["QryEndPos"], row["RefStartPos"], row["RefEndPos"], row["Orientation"],
                               row["Confidence"], row["QryLen"], row["Alignment"])
