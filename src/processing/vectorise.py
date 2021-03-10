from typing import Iterable, List
from itertools import accumulate, zip_longest


def vectoriseSegments(segments: Iterable[int], resolution=100):
    positions = list(accumulate(segments, lambda sumOfSegments, segment: sumOfSegments + segment))
    return vectorise(positions, resolution)


def vectorise(positions: Iterable[int], resolution=100):
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


def blur(vector: List[int], radius: int):
    shifts = range(1, radius + 1)
    shiftedVectors = [vector]
    for shift in shifts:
        shiftedVectors.append(vector[shift:])
        shiftedVectors.append(shift * [0] + vector)

    return [1 if any(position) else 0 for position in zip_longest(*shiftedVectors, fillvalue=0)][0: len(vector)]
