
from validator import Validator


def workerFunction(input):
    (alignment, reference, query, resolution) = input
    result = query.correlate(reference, reverseStrand=alignment.reverseStrand)
    validator = Validator(resolution)
    isValid = validator.validate(result, alignment)

    return 1 if isValid else 0
