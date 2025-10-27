from typing import Iterable

from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segments import AlignmentSegment


class MockSegmentChainer(SegmentChainer):
    def __init__(self):
        super().__init__(None)  # type: ignore

    def chain(self, segments: Iterable[AlignmentSegment], qPositions):
        return list(segments)
