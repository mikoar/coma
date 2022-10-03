import random
from typing import List

import pytest

from src.correlation.bionano_alignment import BionanoAlignment
from src.diagnostic.alignment_comparer import AlignmentRowComparer
from src.diagnostic.xmap_alignment import XmapAlignedPair, XmapAlignmentPosition, XmapAlignedPairWithDistance
from tests.test_doubles.alignment_segment_stub import ScoredAlignedPairStub


def __pair(referenceSiteId: int, querySiteId: int):
    return XmapAlignedPairWithDistance(XmapAlignmentPosition(referenceSiteId, 0), XmapAlignmentPosition(querySiteId, 0),
                                       random.randint(0, 1000))


@pytest.mark.parametrize("referencePairs, queryPairs, alignment1Coverage, alignment2Coverage", [
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     [__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     1., 1.),
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     [],
     0., 1.),
    ([], [__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     1., 0.),
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     [__pair(2, 1), __pair(3, 2), __pair(4, 3)],
     0., 0.),
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3), __pair(4, 4)],
     [__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     0.75, 1.),
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     [__pair(1, 1), __pair(2, 2), __pair(3, 3),
      __pair(4, 4)],
     1., 0.75),
], ids=["identical", "empty query2", "empty query1", "shifted", "one missing in query2", "one missing in query1"])
def test_queryCoverage(referencePairs: List[XmapAlignedPair], queryPairs: List[XmapAlignedPair],
                       alignment1Coverage, alignment2Coverage):
    referenceAlignment = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, referencePairs)
    actualAlignment = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, queryPairs)

    result = AlignmentRowComparer().compare(referenceAlignment, actualAlignment)

    assert alignment1Coverage == result.alignment1Coverage
    assert alignment2Coverage == result.alignment2Coverage


@pytest.mark.parametrize("referencePairs, queryPairs, identity", [
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     [__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     1.),
    ([__pair(1, 1), __pair(2, 2), __pair(3, 3)],
     [__pair(1, 2), __pair(2, 3), __pair(3, 4)],
     0.),
    ([__pair(1, 1), __pair(2, 3)],
     [__pair(1, 1), __pair(3, 2)],
     0.5),

], ids=["identical", "different", "half identical "])
def test_identity(referencePairs: List[XmapAlignedPair], queryPairs: List[ScoredAlignedPairStub], identity):
    referenceAlignment = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, referencePairs)
    actualAlignment = BionanoAlignment(1, 1, 1, 0, 99, 0, 99, False, 123, "1M", 100, 100, queryPairs)

    result = AlignmentRowComparer().compare(referenceAlignment, actualAlignment)

    assert identity == result.identity


if __name__ == '__main__':
    pytest.main(args=[__file__])