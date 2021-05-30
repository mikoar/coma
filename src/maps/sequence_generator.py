
from processing.vectorise import blur, vectorise, vectoriseSegments


class SequenceGenerator:
    def __init__(self, resolution: int, blurRadius: int) -> None:
        self.resolution = resolution
        self.blurRadius = blurRadius

    def segmentsToSequence(self, segments):
        vector = list(vectoriseSegments(segments, self.resolution))
        return blur(vector, self.blurRadius)

    def positionsToSequence(self, positions):
        vector = list(vectorise(positions, self.resolution))
        return blur(vector, self.blurRadius)
