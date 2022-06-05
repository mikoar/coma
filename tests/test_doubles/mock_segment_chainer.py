from typing import Iterable

from src.alignment.segment_chainer import SegmentChainer
from src.alignment.segments import AlignmentSegment


class MockSegmentChainer(SegmentChainer):
    def chain(self, segments: Iterable[AlignmentSegment]):
        return list(segments)
