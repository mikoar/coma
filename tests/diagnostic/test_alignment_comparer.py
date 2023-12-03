from typing import Tuple
from unittest.mock import Mock, call

import pytest

from src.diagnostic.alignment_comparer import AlignmentRowComparer, AlignmentRowComparison, AlignmentComparer, \
    AlignmentRowComparisonResultType
from src.diagnostic.benchmark_alignment import BenchmarkAlignment


class __BenchmarkAlignmentStub(BenchmarkAlignment):
    def __init__(self, queryId: int, referenceId: int):
        self.queryId = queryId
        self.referenceId = referenceId
        self.alignedPairs = []


def __getSut() -> Tuple[AlignmentComparer, Mock]:
    rowComparer: AlignmentRowComparer = Mock(spec=AlignmentRowComparer)
    rowCompareMock = Mock(return_value=AlignmentRowComparison(
        AlignmentRowComparisonResultType.BOTH, BenchmarkAlignment.null, BenchmarkAlignment.null, [], [], 0., 0., 0.))
    rowComparer.compare = rowCompareMock
    return AlignmentComparer(rowComparer), rowCompareMock


def test_matchesAlignmentsByQueryAndReferenceIds():
    alignments1 = [__BenchmarkAlignmentStub(1, 1), __BenchmarkAlignmentStub(2, 1), __BenchmarkAlignmentStub(3, 1)]
    alignments2 = [__BenchmarkAlignmentStub(2, 1), __BenchmarkAlignmentStub(3, 1), __BenchmarkAlignmentStub(4, 1)]

    sut, rowComparerMock = __getSut()
    result = sut.compare(alignments1, alignments2)

    rowComparerMock.assert_has_calls([
        call(alignments1[1], alignments2[0]),
        call(alignments1[2], alignments2[1])
    ])

    assert AlignmentRowComparison.alignment1Only(alignments1[0]) in result.rows
    assert AlignmentRowComparison.alignment2Only(alignments2[2]) in result.rows


if __name__ == '__main__':
    pytest.main(args=[__file__])
