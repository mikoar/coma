import itertools
import pandas
import re
from typing import List, Callable
from functools import cache
from pandas.core.frame import DataFrame
from pandas.core.series import Series
from sequence_generator import SequenceGenerator
from optical_map import OpticalMap, ReferenceOpticalMap


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

    def readQueries(self, filePath: str,  moleculeIds: List[int] = []) -> List[OpticalMap]:
        return self.__read(OpticalMap, filePath, moleculeIds)

    def readQuery(self, filePath: str,  moleculeId: int) -> OpticalMap:
        return self.__read(OpticalMap, filePath, [moleculeId])[0]

    @cache
    def readReference(self, filePath: str, chromosome: int = 1):
        return self.__read(ReferenceOpticalMap, filePath, [chromosome])[0]

    def __read(self, mapConstructor: Callable, filePath, moleculeIds=None):
        maps = self.reader.readFile(filePath, ["CMapId", "Position"])

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup, self.sequenceGenerator, mapConstructor)
        return opticalMaps.tolist()

    @staticmethod
    def __parseCmapRowsGroup(group,  sequenceGenerator: SequenceGenerator, mapConstructor: Callable):
        moleculeId = group["CMapId"].iloc[0]
        positions = group["Position"].sort_values().tolist()
        sequence = sequenceGenerator.positionsToSequence(positions)
        return mapConstructor(moleculeId, sequence, positions, sequenceGenerator.resolution)


class Alignment:
    def __init__(self, queryId, refId, refStart, refEnd, orientation, confidence) -> None:
        self.queryId = queryId
        self.chromosome = refId
        self.refStartPosition = int(refStart)
        self.refEndPosition = int(refEnd)
        self.reverseStrand = orientation == "-"
        self.confidence = confidence


class AlignmentReader:
    def __init__(self) -> None:
        self.reader = BionanoFileReader()

    def readAlignments(self, filePath: str) -> List[Alignment]:
        alignments = self.reader.readFile(filePath,
                                          ["QryContigID", "RefContigID", "RefStartPos",
                                           "RefEndPos", "Orientation", "Confidence"])
        return alignments.apply(self.__parseRow, axis=1).tolist()

    @staticmethod
    def __parseRow(row: Series):
        return Alignment(row["QryContigID"], row["RefContigID"], row["RefStartPos"],
                         row["RefEndPos"], row["Orientation"], row["Confidence"])
