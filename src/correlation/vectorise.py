from itertools import accumulate, zip_longest
from typing import Iterable, List

import numpy as np


def vectoriseSegments(segments: Iterable[int], resolution=100):
    positions = np.array(np.fromiter(accumulate(segments, lambda sumOfSegments, segment: sumOfSegments + segment), dtype=int))
    return vectorisePositions(positions, resolution)


def vectorisePositions(positions: Iterable[int], resolution: int = 100):
    if not isinstance(resolution, int) or resolution < 1:
        raise ValueError(resolution)

    window_start = 0
    window_end = window_start + resolution
    for position in positions:
        if position < window_start:
            continue
        while position >= window_end:
            window_start += resolution
            window_end += resolution
            yield 0

        yield 1
        window_start += resolution
        window_end += resolution


def blur(vector: List[int], radius: int) -> np.ndarray:
    if not isinstance(radius, int) or radius < 0:
        raise ValueError(radius)

    shifts = range(1, radius + 1)
    shiftedVectors = [vector]
    for shift in shifts:
        shiftedVectors.append(vector[shift:])
        shiftedVectors.append(shift * [0] + vector)

    return np.array([1 if any(position) else 0 for position in zip_longest(*shiftedVectors, fillvalue=0)][0: len(vector)])
