from typing import List, TextIO, Iterable

from src.correlation.optical_map import OpticalMap
from src.parsers.bionano_file_reader import BionanoFileReader


class CmapReader:
    def __init__(self, reader: BionanoFileReader = None) -> None:
        self.reader = reader or BionanoFileReader()

    def readQueries(self, file: TextIO, moleculeIds: Iterable[int] = None):
        return self.__read(file, moleculeIds or [])

    def readQuery(self, file: TextIO, moleculeId: int):
        return self.__read(file, [moleculeId])[0]

    def readReference(self, file: TextIO, chromosome: int = 1):
        return self.__read(file, [chromosome])[0]

    def readReferences(self, file: TextIO, chromosomes: Iterable[int] = None):
        return self.__read(file, chromosomes or [])

    def __read(self, file: TextIO, moleculeIds: Iterable[int] = None) -> List[OpticalMap]:
        maps = self.reader.readFile(file, ["CMapId", "Position", "LabelChannel"])

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup)
        return [] if opticalMaps.empty else opticalMaps[opticalMaps.notnull()].tolist()

    @staticmethod
    def __parseCmapRowsGroup(group):
        moleculeId = group["CMapId"].iloc[0]
        labelSites = group[group["LabelChannel"] != 0]
        moleculeEndMarker = group[group["LabelChannel"] == 0].iloc[0]
        length = int(moleculeEndMarker["Position"])
        positions = labelSites["Position"].sort_values().tolist()
        return OpticalMap(moleculeId, length, positions) if positions else None
