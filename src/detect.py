from typing import List


def detect(positions: List[int], resolution=1000):
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

# TODO: wstawiać 1 w sąsiednich oknach, zmniejszyć rozdzielczość,
