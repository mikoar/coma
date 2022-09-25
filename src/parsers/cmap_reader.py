from typing import List, TextIO, Iterable

from src.correlation.optical_map import OpticalMap
from src.parsers.bionano_file_reader import BionanoFileReader


class CmapReader:
    def __init__(self) -> None:
        self.reader = BionanoFileReader()

    def readQueries(self, file: TextIO, moleculeIds: Iterable[int] = None):
        return self.__read(file, moleculeIds or [])

    def readQuery(self, file: TextIO, moleculeId: int):
        return self.__read(file, [moleculeId])[0]

    def readReference(self, file: TextIO, chromosome: int = 1):
        return self.__read(file, [chromosome])[0]

    def readReferences(self, file: TextIO, chromosomes: Iterable[int] = None):
        return self.__read(file, chromosomes or [])

    def __read(self, file: TextIO, moleculeIds: Iterable[int] = None) -> List[OpticalMap]:
        maps = self.reader.readFile(file, ["CMapId", "Position", "ContigLength", "LabelChannel"])

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
