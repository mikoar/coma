
from optical_map import Peaks
from validator import Validator


def workerFunction(input):
    (alignment, reference, query, resolution) = input
    result = query.correlate(reference, reverseStrand=alignment.reverseStrand)
    validator = Validator(resolution)
    peaks = Peaks(result)
    isMaxPeakValid = validator.validate(peaks.max, alignment)

    return (1 if isMaxPeakValid else 0, peaks.getRelativeScore(alignment, validator))
