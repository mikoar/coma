import itertools
import pandas
import re
from typing import List, Callable

from pandas.core.frame import DataFrame
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

    def readQueries(self, filePath: str,  moleculeIds: List[int] = []):
        return self.__read(OpticalMap, filePath, moleculeIds)

    def readReference(self, filePath: str, chromosome: int = 1):
        return self.__read(ReferenceOpticalMap, filePath, [chromosome])[0]

    def __read(self, mapConstructor: Callable, filePath, moleculeIds):
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


# class AlignmentReader:
#     def __init__(self) -> None:
#         self.reader = BionanoFileReader()

#     def readReferenceAlignment(self, filePath: str):
#         data = self.reader.readFile(filePath, ["QryContigID", "RefStartPos", "RefEndPos", "Orientation"])
