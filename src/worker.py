
from optical_map import Peaks
from validator import Validator


def workerFunction(input):
    (alignment, reference, query, resolution) = input
    result = query.correlate(reference, reverseStrand=alignment.reverseStrand)
    validator = Validator(resolution)
    peaks = Peaks(result)
    isValid = validator.validate(peaks, alignment)

    return (1 if isValid else 0, peaks.score)
