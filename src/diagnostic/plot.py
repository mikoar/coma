import math
from typing import List, Tuple, Union

import matplotlib.patches as patches
import numpy as np
import seaborn as sns
from matplotlib import cycler, pyplot, rcParams  # type: ignore
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter
from scipy.interpolate import interp1d

from src.correlation.optical_map import CorrelationResult
from src.correlation.peak import Peak

rcParams["lines.linewidth"] = 1
rcParams['axes.prop_cycle'] = cycler(color=["#e74c3c"])


def plotCorrelation(correlationResult: CorrelationResult,
                    expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None) -> Figure:
    plot, _ = __plotCorrelation(correlationResult, correlationResult.resolution, expectedReferenceRanges)
    return plot


def plotRefinedCorrelation(initialCorrelationResult: CorrelationResult,
                           refinedCorrelationResult: CorrelationResult) -> Figure:
    fig, ax = __plotCorrelation(refinedCorrelationResult, refinedCorrelationResult.resolution)
    ax2 = ax.twinx()
    margin = 30000
    leftMargin = min(margin, refinedCorrelationResult.correlationStart)
    rightMargin = min(margin, initialCorrelationResult.correlationEnd - refinedCorrelationResult.correlationEnd)
    initialCorrelationStart = round(
        (refinedCorrelationResult.correlationStart - leftMargin) / initialCorrelationResult.resolution)
    initialCorrelationEnd = round(
        (refinedCorrelationResult.correlationEnd + rightMargin) / initialCorrelationResult.resolution)
    initialCorrelationFragment = initialCorrelationResult.correlation[initialCorrelationStart:initialCorrelationEnd]
    interpolated = interp1d(np.arange(initialCorrelationFragment.size), initialCorrelationFragment, kind='nearest')
    marginInRefinedCorrelationCoordinates = math.ceil((rightMargin + leftMargin) / refinedCorrelationResult.resolution)
    stretched = interpolated(np.linspace(
        0, initialCorrelationFragment.size - 1,
           refinedCorrelationResult.correlation.size + marginInRefinedCorrelationCoordinates))
    xStart = refinedCorrelationResult.correlationStart - leftMargin
    xEnd = refinedCorrelationResult.correlationEnd + rightMargin
    ax2.set_xlim(left=xStart, right=xEnd)
    ax2.set_ylim(bottom=0, top=stretched.max() * 1.1)
    x = range(xStart, xEnd, refinedCorrelationResult.resolution)
    ax2.fill_between(x, stretched, alpha=0.25)
    __plotPeaks(initialCorrelationResult, ax2, maxAnnotations=0, marker="*")
    return fig


def __plotCorrelation(correlationResult: CorrelationResult,
                      resolution: int,
                      expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None) \
        -> Tuple[Figure, Axes]:
    fig: Figure = pyplot.figure(figsize=(20, 4))
    ax: Axes = fig.add_axes([0, 0, 1, 1])
    ax.ticklabel_format(style='plain')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, p: format(int(value), ',')))

    ax.set_xlim(correlationResult.correlationStart, correlationResult.correlationEnd)
    ax.set_ylim(bottom=0, top=correlationResult.correlation.max() * 1.1)

    if expectedReferenceRanges:
        if isinstance(expectedReferenceRanges, tuple):
            expectedReferenceRanges = [expectedReferenceRanges]
        for expectedRange in expectedReferenceRanges:
            __addExpectedStartStopRect(ax, expectedRange, correlationResult)

    x = range(correlationResult.correlationStart, correlationResult.correlationEnd, resolution)
    ax.plot(x, correlationResult.correlation)

    __plotPeaks(correlationResult, ax)
    __markPeakBaseLevel(ax, correlationResult)

    return fig, ax


def __plotPeaks(correlationResult: CorrelationResult, ax: Axes, maxAnnotations=5, marker: str = "x"):
    sortedPeaks = sorted([peak for peak in correlationResult.peaks if peak.height >= correlationResult.peakBaseLevel],
                         key=lambda p: p.score, reverse=True)
    if not sortedPeaks:
        return

    maxPeak = sortedPeaks[0]
    ax.plot(maxPeak.position, maxPeak.height, marker, markersize=24, markeredgewidth=4)
    __plotSuboptimalPeaks(ax, sortedPeaks[1:maxAnnotations], 0.6, marker)
    __plotSuboptimalPeaks(ax, sortedPeaks[maxAnnotations:], 0.3, marker)

    peaksToAnnotate = sortedPeaks[:maxAnnotations]
    for i, peak in enumerate(peaksToAnnotate):
        ax.annotate(f"({int(peak.position / 1000):,}K, {peak.height:.2f}), s:{peak.score:.3f}",
                    (peak.position, peak.height), rotation=-45, ha="center", va="top")


def __plotSuboptimalPeaks(ax, peaks: List[Peak], alpha: float, marker: str):
    if peaks:
        ax.plot([p.position for p in peaks], [p.height for p in peaks], marker,
                markersize=16, markeredgewidth=4, alpha=alpha)


def __markPeakBaseLevel(ax, correlationResult):
    ax.hlines(correlationResult.peakBaseLevel,
              correlationResult.correlationStart,
              correlationResult.correlationEnd,
              linestyles="--",
              colors="black")


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
