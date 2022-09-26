import pytest

from src.correlation.bionano_alignment import BionanoAlignment
from src.correlation.peak import Peak
from src.diagnostic.validator import Validator


@pytest.mark.parametrize('peakPosition,expected',
                         [(0, True), (10, True), (-10, True), (11, False), (-11, False)])
def test_validator_simple(peakPosition: int, expected: bool):
    qryStart = refStart = 0
    qryEnd = refEnd = 100
    alignment = getBionanoAlignment(qryStart, qryEnd, refStart, refEnd)
    peak = Peak(peakPosition, 1.)
    validator = Validator(1, 10)

    valid = validator.validate(peak, alignment)

    assert valid == expected


@pytest.mark.parametrize('peakPosition,expected',
                         [(0, True), (20, True), (-20, True), (21, False), (-21, False)])
def test_validator_lengthDifference(peakPosition: int, expected: bool):
    qryStart = refStart = 0
    qryEnd = 80
    refEnd = 100
    alignment = getBionanoAlignment(qryStart, qryEnd, refStart, refEnd)
    peak = Peak(peakPosition, 1.)
    validator = Validator(1, 10)

    valid = validator.validate(peak, alignment)

    assert valid == expected


def getBionanoAlignment(qryStart, qryEnd, refStart, refEnd):
    return BionanoAlignment(1, 1, 1, qryStart, qryEnd, refStart, refEnd, '+', 0, '', 0, 0, [])


if __name__ == '__main__':
    pytest.main(args=[__file__])
