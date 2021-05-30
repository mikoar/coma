import itertools
import pandas
import re
from typing import List

from pandas.core.frame import DataFrame
from maps.sequence_generator import SequenceGenerator
from maps.optical_map import OpticalMap, ReferenceOpticalMap


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

    def readQueryMaps(self, filePath: str,  moleculeIds: List[int] = []):
        maps = self.reader.readFile(filePath, ["CMapId", "Position"])

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup, self.sequenceGenerator)
        return opticalMaps.tolist()

    def readReference(self, filePath: str):
        maps = self.reader.readFile(filePath, ["CMapId", "Position"])
        positions = maps["Position"].sort_values().tolist()
        sequence = self.sequenceGenerator.positionsToSequence(positions)
        return ReferenceOpticalMap(sequence, positions)

    @staticmethod
    def __parseCmapRowsGroup(group,  sequenceGenerator: SequenceGenerator):
        moleculeId = group["CMapId"].iloc[0]
        positions = group["Position"].sort_values().tolist()
        sequence = sequenceGenerator.positionsToSequence(positions)
        return OpticalMap(moleculeId, sequence, positions)


class AlignmentReader:
    def __init__(self) -> None:
        self.reader = BionanoFileReader()

    def readReferenceAlignment(self, filePath: str):
        data = self.reader.readFile(filePath, ["QryContigID", "RefStartPos", "RefEndPos", "Orientation"])
