import random
from typing import List, Tuple

import pytest

from src.correlation.bionano_alignment import BionanoAlignment
from src.diagnostic.alignment_comparer import AlignmentRowComparer
from src.diagnostic.benchmark_alignment import BenchmarkAlignmentPosition, \
    BenchmarkAlignedPairWithDistance


def __pair(referenceSiteId: int, querySiteId: int):
    return BenchmarkAlignedPairWithDistance(BenchmarkAlignmentPosition(referenceSiteId, 0), BenchmarkAlignmentPosition(querySiteId, 0),
                                            random.randint(0, 1000))


def __pairs(alignmentPairs: List[Tuple[int, int]]):
    return [__pair(p[0], p[1]) for p in alignmentPairs]


def __compare(alignment1Pairs, alignment2Pairs, combineMultipleQuerySources=False):
    alignment1 = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, __pairs(alignment1Pairs))
    alignment2 = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, __pairs(alignment2Pairs))
    result = AlignmentRowComparer(combineMultipleQuerySources).compare(alignment1, alignment2)
    return result


@pytest.mark.parametrize("alignment1Pairs, alignment2Pairs, alignment1ExclusivePairs, alignment2ExclusivePairs, alignment1Coverage, alignment2Coverage", [
    ([(1, 1), (2, 2), (3, 3)],
     [(1, 1), (2, 2), (3, 3)],
     [], [], 1., 1.),
    ([(1, 1), (2, 2), (3, 3)],
     [],
     [(1, 1), (2, 2), (3, 3)],
     [],
     0., 1.),
    ([], [(1, 1), (2, 2), (3, 3)],
     [], [(1, 1), (2, 2), (3, 3)],
     1., 0.),
    ([(1, 1), (2, 2), (3, 3)],
     [(2, 1), (3, 2), (4, 3)],
     [(1, 1), (2, 2), (3, 3)],
     [(2, 1), (3, 2), (4, 3)],
     0., 0.),
    ([(1, 1), (2, 2), (3, 3), (4, 4)],
     [(1, 1), (2, 2), (3, 3)],
     [(4, 4)],
     [],
     0.75, 1.),
    ([(1, 1), (2, 2), (3, 3)],
     [(1, 1), (2, 2), (3, 3), (4, 4)],
     [],
     [(4, 4)],
     1., 0.75),
    ([(1, 1), (2, 2), (3, 3), (5, 4)],
     [(1, 1), (2, 2), (3, 3), (4, 4)],
     [(5, 4)],
     [(4, 4)],
     0.75, 0.75),
], ids=["identical", "empty alignment2", "empty alignment1", "shifted", "extra pair in alignment1", "one missing in alignment1",
        "differently mapped last query position"])
def test_queryCoverage(
        alignment1Pairs: List[Tuple[int, int]],
        alignment2Pairs: List[Tuple[int, int]],
        alignment1ExclusivePairs: List[Tuple[int, int]],
        alignment2ExclusivePairs: List[Tuple[int, int]],
        alignment1Coverage: float,
        alignment2Coverage: float):
    result = __compare(alignment1Pairs, alignment2Pairs)

    assert result.alignment1Coverage == alignment1Coverage
    assert result.alignment2Coverage == alignment2Coverage
    assert result.alignment1ExclusivePairs == __pairs(alignment1ExclusivePairs)
    assert result.alignment2ExclusivePairs == __pairs(alignment2ExclusivePairs)


@pytest.mark.parametrize("alignment1Pairs, alignment2Pairs, identity", [
    ([(1, 1), (2, 2), (3, 3)],
     [(1, 1), (2, 2), (3, 3)],
     1.),
    ([(1, 1), (2, 2), (3, 3)],
     [(1, 2), (2, 3), (3, 4)],
     0.),
    ([(1, 1), (2, 3)],
     [(1, 1), (3, 2)],
     0.5),

], ids=["identical", "different", "half identical"])
def test_identity(alignment1Pairs: List[Tuple[int, int]], alignment2Pairs: List[Tuple[int, int]], identity):
    result = __compare(alignment1Pairs, alignment2Pairs)
    assert result.identity == identity


def test_combineMultipleQuerySources():
    result = __compare([(1, 1), (4, 2)],
                       [(1, 1), (2, 1), (3, 1), (4, 2)],
                       True)

    assert not result.alignment1ExclusivePairs
    assert not result.alignment2ExclusivePairs
    assert result.identity == 1.


def test_doesNotCombineMultipleQuerySources():
    result = __compare([(1, 1), (4, 2)],
                       [(1, 1), (2, 1), (3, 1), (4, 2)],
                       False)

    assert not result.alignment1ExclusivePairs
    assert result.alignment2ExclusivePairs == __pairs([(2, 1), (3, 1)])
    assert result.identity == 2./3.


if __name__ == '__main__':
    pytest.main(args=[__file__])
