from typing import Tuple
from unittest.mock import Mock, call

import pytest

from src.alignment.alignment_comparer import AlignmentRowComparer, AlignmentRowComparison, AlignmentComparer
from src.correlation.xmap_alignment import XmapAlignment


class __XmapAlignmentStub(XmapAlignment):
    def __init__(self, queryId: int, referenceId: int):
        self.queryId = queryId
        self.referenceId = referenceId
        self.alignedPairs = []


def __getSut() -> Tuple[AlignmentComparer, Mock]:
    rowComparer: AlignmentRowComparer = Mock(spec=AlignmentRowComparer)
    rowCompareMock = Mock(return_value=AlignmentRowComparison(0, 0, [], 0., [], 0., 0.))
    rowComparer.compare = rowCompareMock
    return AlignmentComparer(rowComparer), rowCompareMock


def test_matchesAlignmentsByQueryAndReferenceIds():
    alignments1 = [__XmapAlignmentStub(1, 1), __XmapAlignmentStub(2, 1), __XmapAlignmentStub(3, 1)]
    alignments2 = [__XmapAlignmentStub(2, 1), __XmapAlignmentStub(3, 1), __XmapAlignmentStub(4, 1)]

    sut, rowComparerMock = __getSut()
    result = sut.compare(alignments1, alignments2)

    rowComparerMock.assert_has_calls([
        call(alignments1[1], alignments2[0]),
        call(alignments1[2], alignments2[1])
    ])

    assert AlignmentRowComparison.alignment1Only(1, 1, []) in result.rows
    assert AlignmentRowComparison.alignment2Only(4, 1, []) in result.rows


if __name__ == '__main__':
    pytest.main(args=[__file__])
