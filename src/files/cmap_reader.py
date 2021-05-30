import itertools
import pandas
import re
from typing import List

from pandas.core.frame import DataFrame
from maps.sequence_generator import SequenceGenerator
from maps.optical_map import OpticalMap, ReferenceOpticalMap


class CmapReader:
    def __init__(self, sequenceGenerator: SequenceGenerator) -> None:
        self.sequenceGenerator = sequenceGenerator

    def readQueryMaps(self, filePath: str,  moleculeIds: List[int] = []):
        maps = self.__readCmap(filePath)

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup, self.sequenceGenerator)
        return opticalMaps.tolist()

    def readReference(self, filePath: str):
        maps = self.__readCmap(filePath)
        positions = maps["Position"].sort_values().tolist()
        sequence = self.sequenceGenerator.positionsToSequence(positions)
        return ReferenceOpticalMap(sequence)

    def __readCmap(self, filePath) -> DataFrame:
        return pandas.read_csv(
            filePath, comment="#", delimiter="\t", names=self.__getCmapColumnNames(filePath), usecols=["CMapId", "Position"])  # type: ignore

    def __getCmapColumnNames(self, filePath):
        with open(filePath) as file:
            gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
            header_line = list(itertools.islice(gen, 1))[0].strip()
            names = re.split(r'\s+', header_line)[1:]
        return names

    @staticmethod
    def __parseCmapRowsGroup(group,  sequenceGenerator: SequenceGenerator):
        moleculeId = group["CMapId"].iloc[0]
        positions = group["Position"].sort_values().tolist()
        sequence = sequenceGenerator.positionsToSequence(positions)
        return OpticalMap(moleculeId, sequence)
