# from typing import List

# from src.alignment.alignment_position_scorer import AlignmentPositionScorer


class ChainScorer:
    def __init__(self, 
                 indelOpenPenalty: int,
                 indelExtPenalty: float,
                 minFirstPassScore: int,
                 minSubsequentPassScore: int,
                 isFirstPass: bool=True):
        self.indelOpenPenalty = indelOpenPenalty
        self.indelExtPenalty = indelExtPenalty
        self.minFirstPassScore = minFirstPassScore
        self.minSubsequentPassScore = minFirstPassScore if minSubsequentPassScore is None else minSubsequentPassScore
        self.ultimatePenalty = -10^6
        self.isFirstPass = isFirstPass

    def setFirstPass(self, isFirstPass):
        self.isFirstPass = isFirstPass

    def minAlignmentScore(self):
        return self.minFirstPassScore if self.isFirstPass else self.minSubsequentPassScore

    def scoreIndel(self, indelLength):
        return -self.indelOpenPenalty - self.indelExtPenalty * indelLength

