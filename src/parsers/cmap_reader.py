from typing import List

from src.correlation.optical_map import OpticalMap
from src.parsers.bionano_file_reader import BionanoFileReader


class CmapReader:
    def __init__(self) -> None:
        self.reader = BionanoFileReader()

    def readQueries(self, filePath: str, moleculeIds: List[int] = []):
        return self.__read(filePath, moleculeIds)

    def readQuery(self, filePath: str, moleculeId: int):
        return self.__read(filePath, [moleculeId])[0]

    def readReference(self, filePath: str, chromosome: int = 1):
        return self.__read(filePath, [chromosome])[0]

    def readReferences(self, filePath: str, chromosomes: List[int] = []):
        return self.__read(filePath, chromosomes)

    def __read(self, filePath, moleculeIds=None) -> List[OpticalMap]:
        maps = self.reader.readFile(filePath, ["CMapId", "Position", "ContigLength", "LabelChannel"])

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup)
        return opticalMaps.tolist()

    @staticmethod
    def __parseCmapRowsGroup(group):
        group = group[group["LabelChannel"] != 0]
        moleculeId = group["CMapId"].iloc[0]
        length = int(group["ContigLength"].iloc[0])
        positions = group["Position"].sort_values().tolist()
        return OpticalMap(moleculeId, length, positions)
