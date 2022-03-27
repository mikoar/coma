from typing import List

import pytest

from src.alignment.alignment_comparer import AlignmentComparer
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segments import AlignmentSegment
from src.correlation.bionano_alignment import BionanoAlignment, RefAlignedPair
from tests.test_doubles.alignment_segment_stub import ScoredAlignedPairStub, AlignmentSegmentStub


@pytest.mark.parametrize("referencePairs, segment, query1Coverage, query2Coverage", [
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     AlignmentSegmentStub([ScoredAlignedPairStub(1, 1), ScoredAlignedPairStub(2, 2), ScoredAlignedPairStub(3, 3)]), 1.,
     1.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)], AlignmentSegmentStub([]), 0., 1.),
    ([], AlignmentSegmentStub([ScoredAlignedPairStub(1, 1), ScoredAlignedPairStub(2, 2), ScoredAlignedPairStub(3, 3)]),
     1., 0.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     AlignmentSegmentStub([ScoredAlignedPairStub(2, 1), ScoredAlignedPairStub(3, 2), ScoredAlignedPairStub(4, 3)]), 0.,
     0.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3), RefAlignedPair(4, 4)],
     AlignmentSegmentStub([ScoredAlignedPairStub(1, 1), ScoredAlignedPairStub(2, 2), ScoredAlignedPairStub(3, 3)]),
     0.75, 1.),
    ([RefAlignedPair(1, 1), RefAlignedPair(2, 2), RefAlignedPair(3, 3)],
     AlignmentSegmentStub([ScoredAlignedPairStub(1, 1), ScoredAlignedPairStub(2, 2), ScoredAlignedPairStub(3, 3),
                           ScoredAlignedPairStub(4, 4)]), 1., 0.75),
], ids=["identical", "empty query2", "empty query1", "shifted", "one missing in query2", "one missing in query1"])
def test_compare(referencePairs: List[RefAlignedPair], segment: AlignmentSegment, query1Coverage, query2Coverage):
    referenceAlignment = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, referencePairs)
    actualAlignment = AlignmentResultRow([segment])

    result = AlignmentComparer().compare(referenceAlignment, actualAlignment)

    assert query1Coverage == result.query1Coverage
    assert query2Coverage == result.query2Coverage


if __name__ == '__main__':
    pytest.main(args=[__file__])
