from math import floor
from typing import List, Tuple, Union

import matplotlib.axes
import matplotlib.patches as patches
import numpy as np
import seaborn as sns
from matplotlib import cycler, pyplot, rcParams  # type: ignore
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import interp1d

from src.correlation.optical_map import CorrelationResult

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def plotCorrelation(correlationResult: CorrelationResult,
                    expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None):
    plot, _ = __plotCorrelation(correlationResult, correlationResult.resolution, expectedReferenceRanges)
    return plot


def plotRefinedCorrelation(initialCorrelationResult: CorrelationResult, refinedCorrelationResult: CorrelationResult):
    fig, ax = __plotCorrelation(refinedCorrelationResult, refinedCorrelationResult.resolution)
    ax2 = ax.twinx()
    start = round(refinedCorrelationResult.correlationStart / initialCorrelationResult.resolution)
    end = round(refinedCorrelationResult.correlationEnd / initialCorrelationResult.resolution)
    correlationFragment = initialCorrelationResult.correlation[start:end]
    interpolated = interp1d(np.arange(correlationFragment.size), correlationFragment, kind='nearest')
    stretched = interpolated(np.linspace(0, correlationFragment.size - 1, refinedCorrelationResult.correlation.size))
    x = range(refinedCorrelationResult.correlationStart, refinedCorrelationResult.correlationEnd,
              refinedCorrelationResult.resolution)
    ax2.fill_between(x, stretched, alpha=0.25)
    return fig


def __plotCorrelation(correlationResult: CorrelationResult, resolution: int,
                      expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None):
    fig = pyplot.figure(figsize=(20, 4))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.ticklabel_format(style='plain')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, p: format(int(value), ',')))

    ax.set_xlim(correlationResult.correlationStart, correlationResult.correlationEnd)

    if expectedReferenceRanges:
        if isinstance(expectedReferenceRanges, tuple):
            expectedReferenceRanges = [expectedReferenceRanges]
        for expectedRange in expectedReferenceRanges:
            __addExpectedStartStopRect(ax, expectedRange, correlationResult)

    x = range(correlationResult.correlationStart, correlationResult.correlationEnd, resolution)
    ax.plot(x, correlationResult.correlation)

    __plotPeaks(correlationResult, ax)
    return fig, ax


def __plotPeaks(peaks: CorrelationResult, ax: matplotlib.axes.Axes):
    maxPeak = peaks.maxPeak
    if not maxPeak:
        return

    peaksExceptMax = [peak for peak in peaks.peaks if peak.position != maxPeak.position]
    ax.plot(maxPeak.position, maxPeak.height, "x", markersize=24, markeredgewidth=4)
    if peaksExceptMax:
        ax.plot([p.position for p in peaksExceptMax], peaks.correlation[[
            floor((p.position - peaks.correlationStart) / peaks.resolution) for p in peaksExceptMax]], "x",
                markersize=16, markeredgewidth=4, alpha=0.5)

    labelVerticalPositionIncrement = maxPeak.height / 16
    for i, peak in enumerate(peaks.peaks):
        ax.annotate(f"({int(peak.position / 1000):,}K, {peak.height:.2f})",
                    (peak.position, labelVerticalPositionIncrement * (i % 2)))


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


def plotHeatMap(arr, fileName, x, y):
    ax = sns.heatmap(arr, linewidth=0.5, annot=True, xticklabels=x, yticklabels=y, fmt='.2f')
    ax.get_figure().savefig(fileName)
