from __future__ import annotations

from dataclasses import dataclass
from itertools import groupby
from math import ceil
from typing import List

import numpy as np
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle

from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.optical_map import OpticalMap, PositionWithSiteId
from src.correlation.peak import Peak


@dataclass
class Options:
    limitQueryToAlignedArea = False
    centerPeaks = True


class AlignmentPlot:
    def __init__(self, reference: OpticalMap, query: OpticalMap, alignment: AlignmentResultRow,
                 options: Options = None):
        self.options = options or Options()
        self.figure: Figure = pyplot.figure(figsize=(24, 16))
        self.axes: Axes = self.figure.add_axes([0, 0, 1, 1])
        self.__setLimits(alignment, query)
        self.__plotReference(alignment, reference)
        self.__plotQuery(alignment, query)
        self.__plotSegments(alignment)
        self.axes.legend()

    def __setLimits(self, alignment: AlignmentResultRow, query: OpticalMap):
        offset = alignment.queryLength / 10
        self.yMin = (alignment.queryStartPosition if self.options.limitQueryToAlignedArea else 0) - offset
        drawPeakPositionsWithOffset = list(
            map(lambda s: self.__getDrawPeakPosition(s.peak, alignment) - offset, alignment.segments))
        self.xMin = min([alignment.referenceStartPosition - offset] + drawPeakPositionsWithOffset)
        self.xMax = alignment.referenceEndPosition + offset
        self.axes.set_xlim(self.xMin - offset, self.xMax)
        self.yMax = alignment.queryEndPosition if self.options.limitQueryToAlignedArea else query.length
        self.axes.set_ylim(self.yMin - offset, self.yMax + offset)
        self.axes.set_aspect("equal")

    def __plotReference(self, alignment: AlignmentResultRow, reference: OpticalMap):
        self.axes.set_xlabel(f"Reference {reference.moleculeId}")
        referenceLabelsInScope = [p for p in reference.getPositionsWithSiteIds() if
                                  self.xMin <= p.position <= alignment.referenceEndPosition]
        self.axes.plot([r.position for r in referenceLabelsInScope],
                       np.repeat(self.yMin, len(referenceLabelsInScope)),
                       marker="|",
                       markersize=16,
                       markeredgewidth="2", linewidth=16, markeredgecolor="black",
                       color="yellow")

        for r in self.__skipDensePositions(referenceLabelsInScope):
            self.axes.annotate(r.siteId, (r.position, self.yMin),
                               textcoords="offset points",
                               xytext=(0, -10),
                               ha="center",
                               va="top",
                               rotation=90)

    def __plotQuery(self, alignment: AlignmentResultRow, query: OpticalMap):
        self.axes.set_ylabel(f"Query {query.moleculeId}")
        queryLabelsInScope = [p for p in query.getPositionsWithSiteIds() if
                              not self.options.limitQueryToAlignedArea
                              or self.yMin <= p.position <= alignment.queryEndPosition]

        self.axes.plot(np.repeat(self.xMin, len(queryLabelsInScope)), [q.position for q in queryLabelsInScope],
                       marker="_",
                       markersize=16,
                       markeredgewidth="2",
                       linewidth=16,
                       markeredgecolor="black",
                       color="lime")

        for q in self.__skipDensePositions(queryLabelsInScope):
            self.axes.annotate(q.siteId, (self.xMin, q.position),
                               textcoords="offset points",
                               xytext=(-15, 0),
                               ha="right",
                               va="center")

    def __plotSegments(self, alignment: AlignmentResultRow):
        groupedSegments = groupby(alignment.segments, lambda s: s.peak)
        colors = self.__getContrastingColors(len(alignment.segments))

        for (peakNumber, (peak, segments)), color in zip(enumerate(groupedSegments), colors):
            self.__drawPeak(alignment, color, peak)
            for segmentNumber, segment in enumerate(segments):
                x = list(map(lambda p: p.reference.position, segment.alignedPositions))
                y = list(map(lambda p: p.query.position, segment.alignedPositions))
                self.__plotSegment(color, peakNumber, segmentNumber, x, y)
                self.__drawGrid(x, y)

    def __drawPeak(self, alignment: AlignmentResultRow, color, peak: Peak):
        drawPeakPosition = self.__getDrawPeakPosition(peak, alignment)
        peakRectangle = Rectangle((drawPeakPosition, self.yMin), peak.width, self.yMax - self.yMin,
                                  color=color,
                                  alpha=0.2)
        self.axes.add_patch(peakRectangle)
        self.axes.annotate(f"Peak height: {peak.height}", (drawPeakPosition, self.yMax), rotation=90, va="top")

    def __plotSegment(self, color, peakNumber: int, segmentNumber: int, x, y):
        self.axes.plot(x, y,
                       label=f" peak {peakNumber + 1} segment {segmentNumber + 1}",
                       marker="+",
                       markersize=8,
                       markeredgecolor=color,
                       linewidth=2,
                       color=color)

    def __drawGrid(self, x, y):
        self.axes.vlines(x, self.yMin, y, linestyles="--", colors="black", alpha=0.5, linewidth=0.5)
        self.axes.hlines(y, self.xMin, x, linestyles="--", colors="black", alpha=0.5, linewidth=0.5)

    def __skipDensePositions(self, positions: List[PositionWithSiteId]):
        minDistance = (self.xMax - self.xMin) / 150
        iterator = iter(positions)
        currentPosition = next(iterator)
        previousPosition: PositionWithSiteId | None = None
        while currentPosition is not None:
            if not previousPosition or currentPosition.position - previousPosition.position >= minDistance:
                yield currentPosition
                previousPosition = currentPosition
            currentPosition = next(iterator, None)

    def __getDrawPeakPosition(self, peak: Peak, alignment: AlignmentResultRow):
        return peak.leftProminenceBasePosition + (alignment.queryLength / 2 if self.options.centerPeaks else 0)

    @staticmethod
    def __getContrastingColors(count: int):
        colorMap = pyplot.get_cmap("plasma")
        contrasting = np.column_stack((np.linspace(0, .5, ceil(count / 2)),
                                       np.linspace(.5, 1, ceil(count / 2)))).flatten()
        return colorMap(contrasting)
