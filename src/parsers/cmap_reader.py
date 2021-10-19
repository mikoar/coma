from typing import List

from src.correlation.optical_map import VectorisedOpticalMap
from src.correlation.sequence_generator import SequenceGenerator
from src.parsers.bionano_file_reader import BionanoFileReader


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
