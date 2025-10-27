from itertools import zip_longest
from typing import List

import numpy as np


# def vectorisePositions(positions: List[int], resolution: int = 100, start: int = 0, end: int = None):
#     if not isinstance(resolution, int) or resolution < 1:
#         raise ValueError(resolution)
#     end = end or positions[-1]
#     window_start = start
#     window_end = window_start + resolution
#     for position in positions:
#         if position < window_start:
#             continue
#         while position >= window_end:
#             window_start += resolution
#             window_end += resolution
#             yield 0
#             if window_start > end:
#                 return

#         yield 1
#         window_start += resolution
#         window_end += resolution

def vectorisePositions(positions: List[int], resolution: int = 100, start: int = 0, end: int = None):
    if not isinstance(resolution, int) or resolution < 1:
        raise ValueError(resolution)
    vectSize = (end-start+1)//resolution
    vectPositions = set((pos-start)//resolution for pos in positions)
    return (1 if i in vectPositions else 0 for i in range(vectSize))


def blur(vector: List[int], radius: int) -> np.ndarray:
    if not isinstance(radius, int) or radius < 0:
        raise ValueError(radius)

    shifts = range(1, radius + 1)
    shiftedVectors = [vector]
    for shift in shifts:
        shiftedVectors.append(vector[shift:])
        shiftedVectors.append(shift * [0] + vector)

    return np.array(
        [1 if any(position) else 0 for position in zip_longest(*shiftedVectors, fillvalue=0)][0: len(vector)])
