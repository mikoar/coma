
from validator import Validator


def doStuff(resolution, input):
    (alignment, reference, query) = input
    result = query.correlate(reference, reverseStrand=alignment.reverseStrand)
    validator = Validator(resolution)
    isValid = validator.validate(result, alignment)

    return 1 if isValid else 0
