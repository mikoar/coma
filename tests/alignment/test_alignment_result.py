from __future__ import annotations

from typing import List, Tuple

import pytest
from mock import Mock, call

from src.alignment.aligned_pair import AlignedPair
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.region_score_penalties import RegionScorePenalties


def test_getRegionScores_appliesPenalties_1pair():
    pairs = [AlignedPair(1, 1)]
    alignmentResult = AlignmentResultRow(pairs)
    penalties = __getPenaltiesMock()
    penalties.getDistancePenalty = Mock(return_value=2)
    penalties.getUnmatchedLabelPenalty = Mock(return_value=5)

    perfectMatchScore = 100
    regionScores = alignmentResult.getRegionScores(penalties, perfectMatchScore)

    assert len(regionScores.scores) == 1
    assert regionScores.scores == [100 - 2 - 5]
    penalties.getDistancePenalty.assert_called_once_with(None, pairs[0])
    penalties.getUnmatchedLabelPenalty.assert_called_once_with(None, pairs[0])


def test_getRegionScores_appliesPenalties_2pairs():
    pairs = [AlignedPair(1, 1), AlignedPair(1, 1)]
    alignmentResult = AlignmentResultRow(pairs)
    penalties = __getPenaltiesMock()
    penalties.getDistancePenalty = Mock(return_value=2)
    penalties.getUnmatchedLabelPenalty = Mock(return_value=5)

    perfectMatchScore = 100
    regionScores = alignmentResult.getRegionScores(penalties, perfectMatchScore)

    assert len(regionScores.scores) == 2
    assert regionScores.scores == [100 - 2 - 5] * 2
    calls = [call(None, pairs[0]), call(pairs[0], pairs[1])]
    penalties.getDistancePenalty.assert_has_calls(calls)
    penalties.getUnmatchedLabelPenalty.assert_has_calls(calls)


@pytest.mark.parametrize("pairs, expected", [
    ([(1, 1), (2, 2), (3, 3)], "3M"),
    ([(1, 9), (2, 8), (3, 7)], "3M"),
    ([(20, 42), (21, 43), (22, 44), (23, 45)], "4M"),
    ([(1, 1), (2, 1)], "1D1M"),
    ([(1, 1), (3, 1)], "2D1M"),
    ([(1, 1), (3, 2)], "1D1M"),
    ([(1, 1), (4, 2)], "2D1M"),
    ([(1, 1), (2, 1), (4, 2), (5, 3)], "1D1M1D1M"),
    ([(1, 1), (2, 3)], "1M1I1M"),
    ([(1, 3), (2, 1)], "1M1I1M"),
    ([(1, 1), (2, 4)], "1M2I1M"),
    ([(22904, 1), (22905, 2), (22906, 3), (22907, 4), (22908, 5), (22909, 6), (22910, 7), (22911, 8), (22912, 9), (22913, 10), (22914, 11), (22915, 12), (22916, 13), (22917, 14), (22918, 15), (22920, 15), (22921, 16), (22922, 17), (22923, 17), (22924, 18), (22925, 19), (22926, 20), (22927, 21), (22928, 22), (22929, 23), (22930, 24), (22931, 24), (22932, 25), (22933, 26), (22934, 27), (22935, 28), (22936, 29), (22937, 30), (22938, 31), (22939, 32), (22940, 33), (22941, 34), (22942, 35), (22943, 36), (22944, 37), (22945, 38), (22946, 38), (22947, 39), (22948, 40), (22949, 41), (22950, 42), (22951, 43), (22952, 44), (22953, 44), (22954, 45), (22955, 46), (22956, 46), (22957, 47), (22958, 48), (22959, 49), (22960, 50), (22961, 51), (22962, 52), (22963, 52), (22964, 53), (22966, 54), (22967, 55), (22969, 56), (22970, 57), (22971, 58), (22973, 59), (22974, 60), (22975, 60), (22976, 61), (22977, 61), (22978, 62), (22980, 63), (22981, 64), (22982, 65), (22983, 66), (22984, 67), (22985, 68), (22986, 69), (22987, 70), (22988, 71), (22989, 72), (22990, 73), (22991, 74), (22992, 74), (22993, 75), (22994, 76), (22995, 76), (22996, 77), (22997, 78), (22998, 79), (22999, 80), (23000, 81), (23001, 82), (23002, 83), (23003, 84), (23004, 85), (23005, 86), (23006, 87), (23007, 88), (23008, 89), (23009, 90), (23010, 91), (23011, 92), (23012, 93), (23013, 94), (23014, 95), (23015, 96), (23016, 97), (23017, 98), (23018, 99), (23019, 100), (23020, 101), (23021, 102), (23022, 103), (23023, 104), (23024, 105), (23025, 106), (23027, 106), (23028, 107), (23029, 108), (23030, 109), (23031, 110), (23032, 111), (23033, 111), (23034, 112), (23035, 112), (23036, 113), (23037, 114), (23038, 115), (23039, 115), (23040, 116), (23041, 117), (23042, 118), (23043, 118), (23044, 120), (23045, 121), (23046, 122), (23047, 123), (23048, 124), (23049, 125), (23050, 126), (23051, 127), (23052, 128), (23053, 129), (23054, 130), (23055, 131), (23056, 132), (23057, 132), (23058, 133), (23060, 133), (23061, 134), (23062, 135), (23063, 136), (23064, 137), (23065, 138), (23066, 139), (23067, 140), (23068, 141), (23069, 142), (23070, 143), (23071, 144), (23072, 145), (23073, 146), (23074, 147), (23075, 148), (23077, 148), (23078, 149), (23079, 150), (23080, 151), (23081, 152), (23082, 153), (23083, 154), (23084, 155), (23085, 156), (23086, 157), (23087, 158), (23088, 159), (23089, 160), (23090, 161), (23091, 162), (23092, 163), (23093, 164), (23094, 165), (23095, 166), (23096, 167), (23097, 168), (23098, 169), (23099, 170), (23100, 171), (23101, 172), (23102, 173), (23103, 174), (23104, 175), (23105, 176), (23106, 177), (23107, 178), (23108, 178), (23109, 179), (23110, 179), (23111, 180), (23112, 181), (23125, 205), (23126, 206), (23146, 221), (23147, 222), (23148, 223), (23149, 224), (23150, 225), (23151, 226), (23152, 227), (23153, 228), (23154, 229), (23155, 230), (23156, 231), (23157, 232), (23158, 233), (23159, 234), (23160, 235), (23161, 236), (23162, 237), (23163, 238), (23164, 239), (23165, 240), (23166, 241), (23167, 241), (23169, 242), (23170, 243), (23171, 244), (23172, 245), (23173, 246), (23174, 247)],
     "14M2D2M1D7M1D14M1D6M1D2M1D6M1D2M1D2M1D3M1D1M1D1M1D2M1D11M1D2M1D30M2D5M1D1M1D3M1D3M1D1M1I12M1D1M2D15M2D30M1D1M1D3M23I12D2M14I19D20M1D1M1D6M")
])
def test_cigarString(pairs: List[Tuple[int, int]], expected: str):
    row = AlignmentResultRow(list(map(lambda p: AlignedPair(p[0], p[1]), pairs)))

    assert row.cigarString == expected


def __getPenaltiesMock() -> RegionScorePenalties | Mock:
    return Mock(spec=RegionScorePenalties)


if __name__ == '__main__':
    pytest.main(args=[__file__])
