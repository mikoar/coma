from typing import List

import pytest

from src.alignment.alignment_comparer import AlignmentComparer
from src.alignment.alignment_position import AlignedPair
from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.bionano_alignment import BionanoAlignment, RefAlignedPair
from tests.alignment.test_alignment_position import TestAlignedPair


@pytest.mark.parametrize("pairs1, pairs2, query1Coverage, query2Coverage", [
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     [TestAlignedPair(1, 1), TestAlignedPair(2, 2), TestAlignedPair(3, 3)], 1., 1.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)], [], 0., 1.),
    ([], [TestAlignedPair(1, 1), TestAlignedPair(2, 2), TestAlignedPair(3, 3)], 1., 0.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     [TestAlignedPair(2, 1), TestAlignedPair(3, 2), TestAlignedPair(4, 3)], 0., 0.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3), RefAlignedPair(4, 4)],
     [TestAlignedPair(1, 1), TestAlignedPair(2, 2), TestAlignedPair(3, 3)], 0.75, 1.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     [TestAlignedPair(1, 1), TestAlignedPair(2, 2), TestAlignedPair(3, 3), TestAlignedPair(4, 4)], 1., 0.75),
], ids=["identical", "empty query2", "empty query1", "shifted", "one missing in query2", "one missing in query1"])
def test_compare(pairs1: List[RefAlignedPair], pairs2: List[AlignedPair], query1Coverage, query2Coverage):
    referenceAlignment = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, pairs1)
    actualAlignment = AlignmentResultRow(pairs2)

    result = AlignmentComparer().compare(referenceAlignment, actualAlignment)

    assert query1Coverage == result.query1Coverage
    assert query2Coverage == result.query2Coverage


if __name__ == '__main__':
    pytest.main(args=[__file__])
