from typing import List, Tuple, Union
import matplotlib.patches as patches
from matplotlib import pyplot
import matplotlib.ticker as ticker
from matplotlib.ticker import FuncFormatter
import numpy as np
from optical_map import CorrelationResult
import seaborn as sns


def __addExpectedStartStopRect(ax, expectedReferenceRange: Tuple[int, int], result: CorrelationResult):
    start = (expectedReferenceRange[0], 0)
    width = expectedReferenceRange[1] - expectedReferenceRange[0]
    height = max(result.correlation)

    rect = patches.Rectangle(start, width, height, edgecolor="none", facecolor="black", alpha=0.2)  # type: ignore
    ax.add_patch(rect)

    ax.text(expectedReferenceRange[0], 0, str(expectedReferenceRange[0]), horizontalalignment='left',
            verticalalignment='top')

    ax.text(expectedReferenceRange[1], 0, str(expectedReferenceRange[1]), horizontalalignment='left',
            verticalalignment='top')


def plotCorrelation(result: CorrelationResult, resolution: int, plotReference=False, expectedReferenceRanges: Union[List[Tuple[int, int]], Tuple[int, int]] = None):
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

    if plotReference:
        maxValue = max(result.correlation)
        ax.plot(result.reference.positions, [maxValue] * len(result.reference.positions), 'bo')
        # ax.fill_between(range(0, len(result.correlation) * resolution, resolution), 0, np.array(result.reference.sequence) * maxValue, alpha=0.2)

    return fig


def __plotPeaks(result, resolution, ax):
    maxPeak = result.peaks.max
    peaksExceptMax = np.array([peak for peak in result.peaks.peaks if peak != maxPeak], dtype=np.int)
    ax.plot(maxPeak * resolution, 1, "x", markersize=24, markeredgewidth=4)
    ax.plot(peaksExceptMax * resolution, result.correlation[peaksExceptMax], "x", markersize=16, markeredgewidth=4, alpha=0.5)


def plotHeatMap(arr, count, x, y):
    plotTitle = f"count_{count}_res_{','.join(str(x) for x in y)}_blur_{','.join(str(x) for x in x)}"

    ax = sns.heatmap(arr, linewidth=0.5, annot=True, xticklabels=x, yticklabels=y)
    ax.get_figure().savefig(
        f"../plots_alignments/heatmap_{plotTitle}.svg")
