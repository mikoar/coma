import itertools
import re
from typing import List

import pandas
from pandas.core.frame import DataFrame
from pandas.core.series import Series

from .alignment import Alignment
from .optical_map import VectorisedOpticalMap
from .sequence_generator import SequenceGenerator


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
        return self.__read(filePath, moleculeIds)

    def readQuery(self, filePath: str,  moleculeId: int):
        return self.__read(filePath, [moleculeId])[0]

    def readReference(self, filePath: str, chromosome: int = 1):
        return self.__read(filePath, [chromosome])[0]

    def readReferences(self, filePath: str, chromosomes: List[int] = []):
        return self.__read(filePath, chromosomes)

    def __read(self, filePath, moleculeIds=None) -> List[VectorisedOpticalMap]:
        maps = self.reader.readFile(filePath, ["CMapId", "Position"])

        if moleculeIds:
            maps = maps[maps["CMapId"].isin(moleculeIds)]

        opticalMaps = maps.groupby("CMapId").apply(self.__parseCmapRowsGroup, self.sequenceGenerator)
        return opticalMaps.tolist()

    @staticmethod
    def __parseCmapRowsGroup(group,  sequenceGenerator: SequenceGenerator):
        moleculeId = group["CMapId"].iloc[0]
        positions = group["Position"].sort_values().tolist()
        sequence = sequenceGenerator.positionsToSequence(positions)
        return VectorisedOpticalMap(moleculeId, sequence, positions, sequenceGenerator.resolution)


class LazyCmapReader(CmapReader):
    def __init__(self, sequenceGenerator: SequenceGenerator) -> None:
        super().__init__(sequenceGenerator)
        self.previousReferences: List[VectorisedOpticalMap] = []
        self.previousQueries: List[VectorisedOpticalMap] = []

    def readReferences(self, filePath: str, chromosomes: List[int] = []):
        newReferenceIds = [c for c in chromosomes if c not in self.__previousReferenceIds()]
        newReferences = super().readReferences(filePath, newReferenceIds)
        referencesRequestedAgain = [r for r in self.previousReferences if r.moleculeId in chromosomes]
        references = newReferences + referencesRequestedAgain
        self.previousReferences = references
        return references

    def readQueries(self, filePath: str,  moleculeIds: List[int] = []):
        newQueryIds = [c for c in moleculeIds if c not in self.__previousQueryIds()]
        newQueries = super().readQueries(filePath, newQueryIds)
        queriesRequestedAgain = [r for r in self.previousQueries if r.moleculeId in moleculeIds]
        queries = newQueries + queriesRequestedAgain
        self.previousQueries = queries
        return queries

    def __previousReferenceIds(self):
        return list(map(lambda r: r.moleculeId, self.previousReferences))

    def __previousQueryIds(self):
        return list(map(lambda r: r.moleculeId, self.previousReferences))


class AlignmentReader:
    def __init__(self) -> None:
        self.reader = BionanoFileReader()

    def readAlignments(self, filePath: str, alignmentIds=None) -> List[Alignment]:
        alignments = self.reader.readFile(filePath,
                                          ["XmapEntryID", "QryContigID", "RefContigID",  "QryStartPos",
                                           "QryEndPos", "RefStartPos",  "RefEndPos", "Orientation", "Confidence", "QryLen"])
        if alignmentIds:
            alignments = alignments[alignments["XmapEntryID"].isin(alignmentIds)]

        return alignments.apply(self.__parseRow, axis=1).tolist()

    @staticmethod
    def __parseRow(row: Series):
        return Alignment(row["XmapEntryID"], row["QryContigID"], row["RefContigID"], row["QryStartPos"],
                         row["QryEndPos"], row["RefStartPos"], row["RefEndPos"], row["Orientation"], row["Confidence"], row["QryLen"])
