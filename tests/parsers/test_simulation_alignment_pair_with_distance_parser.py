import pytest

from src.correlation.optical_map import OpticalMap
from src.diagnostic.benchmark_alignment import BenchmarkAlignedPairWithDistance, BenchmarkAlignmentPosition
from src.parsers.simulation_alignment_pair_parser import SimulationAlignmentPairWithDistanceParser


def test_createsAlignedPairsWithDistance_fromSimulationData():
    reference = OpticalMap(1, 100, [1000, 1300, 1320, 2000])
    query = OpticalMap(10, 100, [100, 200, 320, 500])
    sut = SimulationAlignmentPairWithDistanceParser([reference], [query])
    pairs = sut.parse("1:1;FP;1:2,1:3,FP", 10, 1, False)
    assert pairs == [
        BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(1, 1000), BenchmarkAlignmentPosition(1, 100), 0),
        BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(2, 1300), BenchmarkAlignmentPosition(3, 320), 20),
        BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(3, 1320), BenchmarkAlignmentPosition(3, 320), 0),
    ]


def test_createsAlignedPairsWithDistance_fromSimulationData_reverse():
    reference = OpticalMap(1, 100, [1000, 1100, 1200, 1600])
    query = OpticalMap(10, 100, [90, 220, 300])
    sut = SimulationAlignmentPairWithDistanceParser([reference], [query])
    pairs = sut.parse("1:3;1:2;1:1", 10, 1, True)
    assert pairs == [
        BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(1, 1000), BenchmarkAlignmentPosition(3, 300), 0),
        BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(2, 1100), BenchmarkAlignmentPosition(2, 220), -20),
        BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(3, 1200), BenchmarkAlignmentPosition(1, 90), 10),
    ]


def test_createsAlignedPairsWithDistance_fromSimulationData_empty():
    reference = OpticalMap(1, 100, [1000, 1100, 1200])
    query = OpticalMap(10, 100, [100, 200, 300])
    sut = SimulationAlignmentPairWithDistanceParser([reference], [query])
    pairs = sut.parse("", 10, 1, False)
    assert len(pairs) == 0


if __name__ == '__main__':
    pytest.main(args=[__file__])
