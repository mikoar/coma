from math import floor, log10
from typing import List, Tuple, Union

import matplotlib.axes
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib import cycler, pyplot, rcParams  # type: ignore
from matplotlib.ticker import FuncFormatter

from src.correlation.optical_map import CorrelationResult

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def __addExpectedStartStopRect(ax, expectedReferenceRange: Tuple[int, int], peaks: CorrelationResult):
    start = (expectedReferenceRange[0], 0)
    width = expectedReferenceRange[1] - expectedReferenceRange[0]
    height = peaks.correlation.max()

    rect = patches.Rectangle(start, width, height, edgecolor="none", facecolor="black", alpha=0.2)  # type: ignore
    ax.add_patch(rect)

    ax.text(expectedReferenceRange[0], 0, str(expectedReferenceRange[0]), horizontalalignment='left',
            verticalalignment='top')

    ax.text(expectedReferenceRange[1], 0, str(expectedReferenceRange[1]), horizontalalignment='left',
            verticalalignment='top')


def plotCorrelation(peaks: CorrelationResult, resolution: int,
                    expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None):
    fig = pyplot.figure(figsize=(40, 5))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.ticklabel_format(style='plain')
    length = len(peaks.correlation) * resolution
    ax.xaxis.set_major_locator(ticker.MultipleLocator(10 ** floor(log10(length))))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, p: format(int(x), ',')))

    ax.set_xlim(peaks.correlationStart, peaks.correlationEnd)

    if expectedReferenceRanges:
        if isinstance(expectedReferenceRanges, tuple):
            expectedReferenceRanges = [expectedReferenceRanges]
        for expectedRange in expectedReferenceRanges:
            __addExpectedStartStopRect(ax, expectedRange, peaks)

    x = range(peaks.correlationStart, peaks.correlationEnd, resolution)
    ax.plot(x, peaks.correlation)

    __plotPeaks(peaks, resolution, ax)

    return fig


def __plotPeaks(peaks: CorrelationResult, resolution, ax: matplotlib.axes.Axes):
    maxPeak = peaks.maxPeak
    if not maxPeak:
        return

    peaksExceptMax = [peak for peak in peaks.peaks if peak.position != maxPeak.position]
    ax.plot(maxPeak.positionInReference, maxPeak.height, "x", markersize=24, markeredgewidth=4)
    if peaksExceptMax:
        ax.plot([p.positionInReference for p in peaksExceptMax], peaks.correlation[[
            floor(p.position - peaks.correlationStart / peaks.resolution) for p in peaksExceptMax]], "x", markersize=16, markeredgewidth=4, alpha=0.5)
    for peak in peaks.peaks:
        ax.annotate(f"({int(peak.positionInReference)}, {peak.height})", (peak.positionInReference, 0))


def plotHeatMap(arr, fileName, x, y):
    ax = sns.heatmap(arr, linewidth=0.5, annot=True, xticklabels=x, yticklabels=y, fmt='.2f')
    ax.get_figure().savefig(fileName)
