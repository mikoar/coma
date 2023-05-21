from src.alignment.alignment_position import ScoredAlignedPair
from tests.test_doubles.aligned_pair_builder import AlignedPairBuilder


class ScoredAlignedPairBuilder(AlignedPairBuilder):
    def __init__(self):
        super().__init__()
        self.score = 500.

    def withScore(self, score):
        self.score = score
        return self

    def build(self):
        return ScoredAlignedPair(AlignedPairBuilder.build(self), self.score)
