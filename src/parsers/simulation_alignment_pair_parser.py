from abc import ABC, abstractmethod
from typing import List

from src.correlation.optical_map import OpticalMap
from src.diagnostic.benchmark_alignment import BenchmarkAlignedPair, BenchmarkAlignmentPosition, \
    BenchmarkAlignedPairWithDistance


class BaseSimulationAlignmentPairParser(ABC):
    @abstractmethod
    def parse(self, simulationDetail: str, queryId: int, referenceId: int, reverseStrand: bool):
        pass


class SimulationAlignmentPairParser(BaseSimulationAlignmentPairParser):
    def parse(self, simulationDetail: str, queryId: int, referenceId: int, reverseStrand: bool):
        pairs = [pair for pairs in
                 [SimulationAlignmentPairParser.__parsePair(x, i + 1) for i, x in enumerate(simulationDetail.split(";")) if x != "FP"]
                 for pair in pairs]
        if reverseStrand:
            pairs.reverse()
        return pairs

    @staticmethod
    def __parsePair(simulationDetailOfPosition: str, querySiteId: int):
        if "," in simulationDetailOfPosition:
            splitPositionStrings = simulationDetailOfPosition.split(",")
            return [BenchmarkAlignedPair.create(s.split(":")[1], str(querySiteId)) for s in splitPositionStrings if s != "FP"]
        else:
            return [BenchmarkAlignedPair.create(simulationDetailOfPosition.split(":")[1], str(querySiteId))]


class SimulationAlignmentPairWithDistanceParser(BaseSimulationAlignmentPairParser):
    def __init__(self, references: List[OpticalMap], queries: List[OpticalMap]):
        self.references = references
        self.queries = queries

    def parse(self, simulationDetail: str, queryId: int, referenceId: int, reverseStrand: bool):
        reference = next(r for r in self.references if r.moleculeId == referenceId)
        query = next(q for q in self.queries if q.moleculeId == queryId)

        pairs = [pair for pairs in
                 [self.__parsePair(x, i + 1, reference, query) for i, x in enumerate(simulationDetail.split(";")) if x != "FP"]
                 for pair in pairs]
        if reverseStrand:
            pairs.reverse()

        pairsWithDistance = map(lambda pair: BenchmarkAlignedPairWithDistance.calculateDistance(pair, pairs[0], reverseStrand), pairs)
        return list(pairsWithDistance) if pairs else []

    def __parsePair(self, simulationDetailOfPosition: str, querySiteId: int, reference: OpticalMap, query: OpticalMap):
        if not simulationDetailOfPosition:
            return []
        if "," in simulationDetailOfPosition:
            splitPositionStrings = simulationDetailOfPosition.split(",")
            return [self.__createPairWithDistance(int(s.split(":")[1]), querySiteId, reference, query) for s in splitPositionStrings if s != "FP"]
        else:
            return [self.__createPairWithDistance(int(simulationDetailOfPosition.split(":")[1]), querySiteId, reference, query)]

    @staticmethod
    def __createPairWithDistance(referenceSiteId: int, querySiteId: int, reference: OpticalMap, query: OpticalMap):
        return BenchmarkAlignedPair(
            BenchmarkAlignmentPosition(referenceSiteId, reference.positions[referenceSiteId - 1]),
            BenchmarkAlignmentPosition(querySiteId, query.positions[querySiteId - 1]))
