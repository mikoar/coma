
import pytest
from src.alignment import Alignment
from src.validator import Validator


@pytest.mark.skip
def test_validator_simple():
    qryStart = refStart = 0
    qryEnd = refEnd = 100
    alignment = Alignment(1, 1, 1, qryStart, qryEnd, refStart, refEnd, '+', 1, 100)
    def peaks(): return None
    peaks.max = 50
    validator = Validator(1)

    valid = validator.validate(peaks, alignment)  # type: ignore

    assert valid


@pytest.mark.skip
def test_validator_invalid():
    qryStart = refStart = 0
    qryEnd = refEnd = 100
    alignment = Alignment(1, 1, 1, qryStart, qryEnd, refStart, refEnd, '+', 1, 100)
    def peaks(): return None
    peaks.max = 101
    validator = Validator(1)

    valid = validator.validate(peaks, alignment)  # type: ignore

    assert not valid


@pytest.mark.skip
def test_validator_short_alignment_without_querys_middle():
    qryStart = refStart = 0
    qryEnd = refEnd = 20
    alignment = Alignment(1, 1, 1, qryStart, qryEnd, refStart, refEnd, '+', 1, 100)
    def peaks(): return None
    peaks.max = 50
    validator = Validator(1)

    valid = validator.validate(peaks, alignment)  # type: ignore

    assert valid


@pytest.mark.skip
def test_validator_short_alignment_without_querys_middle_query_at_query_end():
    qryStart = refStart = 80
    qryEnd = refEnd = 100
    alignment = Alignment(1, 1, 1, qryStart, qryEnd, refStart, refEnd, '+', 1, 100)
    def peaks(): return None
    peaks.max = 50
    validator = Validator(1)

    valid = validator.validate(peaks, alignment)  # type: ignore

    assert valid


if __name__ == '__main__':
    pytest.main(args=[__file__])
