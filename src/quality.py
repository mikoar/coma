from optical_map import CorrelationResult
from scipy.signal import find_peaks


class Quality:
    def __init__(self, correlationResult: CorrelationResult) -> None:
        self.correlationResult = correlationResult
        self.peaks, self.peakProperties = find_peaks(correlationResult.correlation,
                                                     height=0.5,
                                                     prominence=0.2,
                                                     distance=(10 ** 8 / correlationResult.query.resolution))
        self.score = self.__getScore()

        self.reverseScore = 1 / self.score if self.score else 1

    def __getScore(self):

        heights = self.peakProperties['peak_heights']
        if len(heights < 2):
            return 0

        twoHighest = sorted(heights, reverse=True)[:2]
        return twoHighest[0] - twoHighest[1]
