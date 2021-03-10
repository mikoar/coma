import matplotlib.patches as patches
from matplotlib import pyplot
import matplotlib.ticker as ticker

from maps.optical_map import CorrelationResult


def __addExpectedStartStopRect(ax, result: CorrelationResult):
    referenceCoordinates = result.query.referenceCoordinates
    start = (referenceCoordinates.start, 0)
    width = referenceCoordinates.stop - referenceCoordinates.start
    height = max(result.correlation)

    rect = patches.Rectangle(start, width, height, edgecolor="none", facecolor="lightgrey")
    ax.add_patch(rect)

    ax.text(referenceCoordinates.start, 0, str(referenceCoordinates.start), horizontalalignment='left',
            verticalalignment='top')

    ax.text(referenceCoordinates.stop, 0, str(referenceCoordinates.stop), horizontalalignment='left',
            verticalalignment='top')


def plotCorrelation(result: CorrelationResult, resolution: int):
    fig = pyplot.figure(figsize=(40, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.ticklabel_format(style='plain')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(100000))

    ax.set_xlim(0, len(result.correlation) * resolution)
    __addExpectedStartStopRect(ax, result)

    ax.plot(range(0, len(result.correlation) * resolution, resolution), result.correlation)
    return fig
