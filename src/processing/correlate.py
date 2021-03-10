
from typing import Iterable
from scipy.signal import correlate


def crossCorrelate(reference: Iterable[int], query: Iterable[int]):
    corr = correlate(reference, query, mode='same')
    return corr
