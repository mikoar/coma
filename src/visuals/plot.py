import matplotlib.patches as patches
from matplotlib import pyplot
import matplotlib.ticker as ticker
import numpy as np

from maps.optical_map import CorrelationResult


def __addExpectedStartStopRect(ax, expectedRange, result: CorrelationResult):
    start = (expectedRange[0], 0)
    width = expectedRange[1] - expectedRange[0]
    height = max(result.correlation)

    rect = patches.Rectangle(start, width, height, edgecolor="none", facecolor="black", alpha=0.2)  # type: ignore
    ax.add_patch(rect)

    ax.text(expectedRange[0], 0, str(expectedRange[0]), horizontalalignment='left',
            verticalalignment='top')

    ax.text(expectedRange[1], 0, str(expectedRange[1]), horizontalalignment='left',
            verticalalignment='top')


def plotCorrelation(result: CorrelationResult, resolution: int, plotReference=False, expectedRange=None):
    fig = pyplot.figure(figsize=(40, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.ticklabel_format(style='plain')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10 ** 7))

    ax.set_xlim(0, len(result.correlation) * resolution)

    if expectedRange:
        __addExpectedStartStopRect(ax, expectedRange, result)

    lenght = len(result.correlation) * resolution
    ax.plot(range(0, lenght, resolution), result.correlation)

    if plotReference:
        maxValue = max(result.correlation)
        ax.plot(result.reference.positions, [maxValue] * len(result.reference.positions), 'bo')
        # ax.fill_between(range(0, len(result.correlation) * resolution, resolution), 0, np.array(result.reference.sequence) * maxValue, alpha=0.2)

    return fig
