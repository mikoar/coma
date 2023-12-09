from abc import ABC, abstractmethod
from typing import List

from src.correlation.optical_map import OpticalMap
from src.diagnostic.benchmark_alignment import BenchmarkAlignedPair, BenchmarkAlignmentPosition, \
    BenchmarkAlignedPairWithDistance


class BaseXmapAlignmentPairParser(ABC):
    @abstractmethod
    def parse(self, alignment: str, queryId: int, referenceId: int, reverseStrand: bool):
        pass


class XmapAlignmentPairParser(BaseXmapAlignmentPairParser):
    def parse(self, alignment: str, queryId: int, referenceId: int, reverseStrand: bool):
        alignmentPairStrings = alignment[:-1].replace('(', '').split(')')
        return list(map(lambda pair: BenchmarkAlignedPair.create(*pair.split(',')), alignmentPairStrings))


class XmapAlignmentPairWithDistanceParser(BaseXmapAlignmentPairParser):
    def __init__(self, references: List[OpticalMap], queries: List[OpticalMap]):
        self.references = references
        self.queries = queries

    def parse(self, alignment: str, queryId: int, referenceId: int, reverseStrand: bool):
        alignmentPairStrings = alignment[:-1].replace('(', '').split(')')
        reference = next(r for r in self.references if r.moleculeId == referenceId)
        query = next(q for q in self.queries if q.moleculeId == queryId)

        def createAlignedPair(pair: str):
            referenceSiteId, querySiteId = map(lambda siteId: int(siteId), pair.split(','))
            referencePosition = BenchmarkAlignmentPosition(referenceSiteId, reference.positions[referenceSiteId - 1])
            queryPosition = BenchmarkAlignmentPosition(querySiteId, query.positions[querySiteId - 1])
            return BenchmarkAlignedPair(referencePosition, queryPosition)

        alignedPairs = list(map(lambda pair: createAlignedPair(pair), alignmentPairStrings))
        alignedPairsWithDistances = map(
            lambda pair: BenchmarkAlignedPairWithDistance.calculateDistance(pair, alignedPairs[0], reverseStrand),
            alignedPairs)
        return list(alignedPairsWithDistances) if alignedPairs else []
