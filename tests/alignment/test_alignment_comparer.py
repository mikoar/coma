from typing import List

import pytest

from src.alignment.aligned_pair import AlignedPair
from src.alignment.alignment_comparer import AlignmentComparer
from src.alignment.alignment_result import AlignmentResult
from src.correlation.alignment import Alignment, RefAlignedPair


@pytest.mark.parametrize("pairs1, pairs2, query1Coverage, query2Coverage", [
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     [AlignedPair(1, 1), AlignedPair(2, 2), AlignedPair(3, 3)], 1., 1.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)], [], 0., 1.),
    ([], [AlignedPair(1, 1), AlignedPair(2, 2), AlignedPair(3, 3)], 1., 0.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     [AlignedPair(2, 1), AlignedPair(3, 2), AlignedPair(4, 3)], 0., 0.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3), RefAlignedPair(4, 4)],
     [AlignedPair(1, 1), AlignedPair(2, 2), AlignedPair(3, 3)], 0.75, 1.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     [AlignedPair(1, 1), AlignedPair(2, 2), AlignedPair(3, 3), AlignedPair(4, 4)], 1., 0.75),
], ids=["identical", "empty query2", "empty query1", "shifted", "one missing in query2", "one missing in query1"])
def test_compare(pairs1: List[RefAlignedPair], pairs2: List[AlignedPair], query1Coverage, query2Coverage):
    referenceAlignment = Alignment(1, 1, 1, 0, 99, 0, 99, False, 123, 100, pairs1)
    actualAlignment = AlignmentResult(0, 99, pairs2)

    result = AlignmentComparer().compare(referenceAlignment, actualAlignment)

    assert query1Coverage == result.query1Coverage
    assert query2Coverage == result.query2Coverage


if __name__ == '__main__':
    pytest.main(args=[__file__])
