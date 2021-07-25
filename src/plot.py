from typing import List, Tuple, Union

import matplotlib.patches as patches
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib import cycler, pyplot, rcParams  # type: ignore
from matplotlib.ticker import FuncFormatter

from .optical_map import CorrelationResult, Peaks

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def __addExpectedStartStopRect(ax, expectedReferenceRange: Tuple[int, int], result: CorrelationResult):
    start = (expectedReferenceRange[0], 0)
    width = expectedReferenceRange[1] - expectedReferenceRange[0]
    height = result.correlation.max()

    rect = patches.Rectangle(start, width, height, edgecolor="none", facecolor="black", alpha=0.2)  # type: ignore
    ax.add_patch(rect)

    ax.text(expectedReferenceRange[0], 0, str(expectedReferenceRange[0]), horizontalalignment='left',
            verticalalignment='top')

    ax.text(expectedReferenceRange[1], 0, str(expectedReferenceRange[1]), horizontalalignment='left',
            verticalalignment='top')


def plotCorrelation(result: CorrelationResult, resolution: int,  expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None):
    fig = pyplot.figure(figsize=(40, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.ticklabel_format(style='plain')
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10 ** 7))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))

    ax.set_xlim(0, len(result.correlation) * resolution)

    if expectedReferenceRanges:
        if isinstance(expectedReferenceRanges, tuple):
            expectedReferenceRanges = [expectedReferenceRanges]
        for expectedRange in expectedReferenceRanges:
            __addExpectedStartStopRect(ax, expectedRange, result)

    lenght = len(result.correlation) * resolution
    x = range(0, lenght, resolution)
    ax.plot(x, result.correlation)

    __plotPeaks(result, resolution, ax)

    return fig


def __plotPeaks(result: CorrelationResult, resolution, ax):
    peaks = Peaks(result)
    maxPeak = peaks.max
    if not maxPeak:
        return

    peaksExceptMax = [peak for peak in peaks.peaks if peak != maxPeak.position]
    ax.plot(maxPeak.positionInReference, 1, "x", markersize=24, markeredgewidth=4)
    ax.plot([p.positionInReference for p in peaksExceptMax], result.correlation[[
            p.position for p in peaksExceptMax]], "x", markersize=16, markeredgewidth=4, alpha=0.5)


def plotHeatMap(arr, fileName, x, y):
    ax = sns.heatmap(arr, linewidth=0.5, annot=True, xticklabels=x, yticklabels=y, fmt='.2f')
    ax.get_figure().savefig(fileName)
