from __future__ import annotations

from dataclasses import dataclass
from itertools import groupby
from math import ceil
from typing import List, Tuple

import numpy as np
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
from matplotlib.ticker import EngFormatter

from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import OpticalMap, PositionWithSiteId, InitialAlignment
from src.correlation.peak import Peak
from src.diagnostic.xmap_alignment import XmapAlignment, XmapAlignedPair


@dataclass
class Options:
    limitQueryToAlignedArea = False
    drawGridForNotAlignedPositions = True


class AlignmentPlot:
    def __init__(self, reference: OpticalMap, query: OpticalMap, alignment: AlignmentResultRow,
                 correlation: InitialAlignment, benchmarkAlignment: XmapAlignment, options: Options = None):
        self.reference = reference
        self.query = query
        self.alignment = alignment
        self.correlation = correlation
        self.benchmarkAlignment = benchmarkAlignment
        self.options = options or Options()
        self.__createFigure()
        self.__setDimensions()
        self.__plotReference()
        self.__plotQuery()
        self.__drawPrimaryPeak()
        self.__plotBenchmarkAlignment()
        self.__plotSegments()
        self.__drawGridForNotAlignedPositions()
        self.axes.legend()

    def __createFigure(self):
        self.figure: Figure = pyplot.figure()
        self.axes: Axes = self.figure.add_axes([0, 0, 1, 1])
        self.axes.set_aspect("equal")
        self.axes.ticklabel_format(style='plain')
        formatter = EngFormatter(unit="bp")
        self.axes.xaxis.set_major_formatter(formatter)
        self.axes.yaxis.set_major_formatter(formatter)

    def __setDimensions(self):
        margin = self.alignment.queryLength / 10
        drawPeakPositions = \
            list(map(lambda s: s.peak.leftProminenceBasePosition, self.alignment.segments)) \
            + [self.correlation.maxPeak.leftProminenceBasePosition]

        self.overlapsWithBenchmark = \
            (not (self.benchmarkAlignment.referenceEndPosition < self.alignment.referenceStartPosition
                  or self.alignment.referenceEndPosition <= self.benchmarkAlignment.referenceStartPosition)) \
                if self.benchmarkAlignment else False

        alignmentStartPositions = [self.alignment.referenceStartPosition]
        alignmentReferenceEndPositions = [self.alignment.referenceEndPosition]

        self.yMinPlot = self.alignment.queryStartPosition if self.options.limitQueryToAlignedArea else 0
        self.yMaxPlot = self.alignment.queryEndPosition if self.options.limitQueryToAlignedArea else self.query.length

        self.yMinAxis = self.yMinPlot - margin
        self.yMinBorder = self.yMinAxis - margin
        self.yMaxAxis = self.yMaxPlot + margin

        if self.overlapsWithBenchmark:
            alignmentStartPositions.append(self.benchmarkAlignment.referenceStartPosition)
            alignmentReferenceEndPositions.append(self.benchmarkAlignment.referenceEndPosition)

        self.xMinPlot = min(alignmentStartPositions + drawPeakPositions)
        self.xMaxPlot = max(alignmentReferenceEndPositions)

        self.xMinAxis = self.xMinPlot - margin
        self.xMinBorder = self.xMinAxis - margin
        self.xMaxAxis = self.xMaxPlot + margin
        self.axes.set_xlim(self.xMinBorder, self.xMaxAxis)
        self.axes.set_ylim(self.yMinBorder, self.yMaxAxis)
        self.axes.set_yticks([y for y in self.axes.get_yticks() if y >= 0])
        self.plotAreaMask = Rectangle((self.xMinAxis, self.yMinAxis),
                                      self.xMaxPlot - self.xMinAxis,
                                      self.yMaxPlot - self.yMinAxis,
                                      facecolor='none', edgecolor='none')

        self.axes.add_patch(self.plotAreaMask)
        self.__setFigureSize()

    def __setFigureSize(self):
        scale = 45_000
        minSize = 10.
        maxSize = 60.
        xSize = (self.xMaxAxis - self.xMinBorder) / scale
        ySize = (self.yMaxAxis - self.yMinBorder) / scale
        xSizeClamped = min(max([xSize, minSize]), maxSize)
        ySizeClamped = min(max([ySize, minSize]), maxSize)
        self.figure.set_size_inches((xSizeClamped, ySizeClamped))

    def __plotReference(self):
        self.axes.set_xlabel(f"Reference {self.reference.moleculeId}")
        refLabelsInScope = [p for p in self.reference.getPositionsWithSiteIds()
                            if self.__isReferencePositionInScope(p.position)]
        self.axes.plot([r.position for r in refLabelsInScope],
                       np.repeat(self.yMinAxis, len(refLabelsInScope)),
                       marker="|",
                       markersize=16,
                       markeredgewidth="2", linewidth=16, markeredgecolor="black",
                       color="yellow")

        for r in self.__skipDensePositions(refLabelsInScope):
            self.axes.annotate(r.siteId, (r.position, self.yMinAxis),
                               textcoords="offset points",
                               xytext=(0, -10),
                               ha="center",
                               va="top",
                               rotation=90)

    def __plotQuery(self):
        self.axes.set_ylabel(f"Query {self.query.moleculeId}")
        queryLabelsInScope = [p for p in self.query.getPositionsWithSiteIds()
                              if not self.options.limitQueryToAlignedArea or self.__isQueryPositionInScope(p.position)]

        self.axes.plot(np.repeat(self.xMinAxis, len(queryLabelsInScope)), [q.position for q in queryLabelsInScope],
                       marker="_",
                       markersize=16,
                       markeredgewidth="2",
                       linewidth=16,
                       markeredgecolor="black",
                       color="lime")

        for q in self.__skipDensePositions(queryLabelsInScope):
            self.axes.annotate(q.siteId, (self.xMinAxis, q.position),
                               textcoords="offset points",
                               xytext=(-15, 0),
                               ha="right",
                               va="center")

    def __plotSegments(self):
        groupedSegments = groupby(self.alignment.segments, lambda s: s.peak)
        colors = self.__getContrastingColors(len(self.alignment.segments))

        for (peakNumber, (peak, segments)), color in zip(enumerate(groupedSegments), colors):
            self.__drawPeak(color, peak)
            for segmentNumber, segment in enumerate(segments):
                x = list(map(lambda p: p.reference.position, segment.alignedPositions))
                y = list(map(lambda p: p.query.position, segment.alignedPositions))
                self.__plotSegment(color, peak, peakNumber, segment, segmentNumber, x, y)
                self.__drawGrid(x, y)

    def __drawPeak(self, color, peak: Peak):
        peakRectangle = Rectangle((peak.leftProminenceBasePosition, self.yMinBorder - 9999999999.),
                                  peak.width,
                                  self.yMaxAxis + 99999999999999.,
                                  color=color,
                                  alpha=0.2,
                                  angle=-45.,
                                  rotation_point=(peak.position, 0))  # type:ignore
        self.axes.add_patch(peakRectangle)
        peakRectangle.set_clip_path(self.plotAreaMask)

    def __drawPrimaryPeak(self):
        peak = self.correlation.maxPeak
        self.axes.plot([peak.position, self.xMaxPlot], [0, self.xMaxPlot - peak.position],
                       linestyle="dashdot",
                       marker=None,
                       color="black")

        peakRectangle = Rectangle((peak.leftProminenceBasePosition, self.yMinBorder - 9999999999.),
                                  peak.width,
                                  self.yMaxAxis + 99999999999999.,
                                  label="primary correlation",
                                  facecolor="wheat",
                                  edgecolor="black",
                                  linestyle="dashdot",
                                  linewidth=0.5,
                                  alpha=0.5,
                                  angle=-45.,
                                  rotation_point=(peak.position, 0))  # type:ignore
        self.axes.add_patch(peakRectangle)
        peakRectangle.set_clip_path(self.plotAreaMask)

    def __plotBenchmarkAlignment(self):
        if not self.benchmarkAlignment:
            return

        x, y = list(zip(*map(lambda p: (p.reference.position, p.query.position), self.benchmarkAlignment.alignedPairs)))
        self.axes.plot(x, y,
                       label=f"benchmark ({len(self.benchmarkAlignment.alignedPairs)} pairs, "
                             f"confidence: {self.benchmarkAlignment.confidence})",
                       color="gray",
                       markeredgecolor="gray",
                       fillstyle="none",
                       marker="o",
                       markersize=8,
                       linewidth=2)

        if self.overlapsWithBenchmark:
            self.__drawGrid(x, y, lineStyle=(0, (1, 5)))

    def __plotSegment(self, color, peak: Peak, peakNumber: int, segment: AlignmentSegment, segmentNumber: int, x, y):
        self.axes.plot(x, y,
                       label=f" peak {peakNumber + 1} (height: {peak.height:.0f}), segment {segmentNumber + 1} "
                             f"({len(segment.alignedPositions)} pairs, score: {segment.segmentScore:.1f})",
                       marker="+",
                       markersize=8,
                       markeredgecolor=color,
                       linewidth=2,
                       color=color)

    def __drawGrid(self, x, y, xMax=None, yMax=None, lineStyle: str | Tuple = "--", label: str = None):
        self.axes.vlines(x, self.yMinAxis, xMax or y, linestyles=lineStyle, colors="gray", linewidth=0.5)
        self.axes.hlines(y, self.xMinAxis, yMax or x, linestyles=lineStyle, colors="gray", linewidth=0.5,
                         label=label)

    def __skipDensePositions(self, positions: List[PositionWithSiteId]):
        minDistance = (self.xMaxPlot - self.xMinPlot) / 150
        iterator = iter(positions)
        currentPosition = next(iterator, None)
        previousPosition: PositionWithSiteId | None = None
        while currentPosition is not None:
            if not previousPosition or currentPosition.position - previousPosition.position >= minDistance:
                yield currentPosition
                previousPosition = currentPosition
            currentPosition = next(iterator, None)

    @staticmethod
    def __getContrastingColors(count: int):
        colorMap = pyplot.get_cmap("plasma")
        contrasting = np.column_stack((np.linspace(0, .5, ceil(count / 2)),
                                       np.linspace(.5, 1, ceil(count / 2)))).flatten()
        return colorMap(contrasting)

    def __drawGridForNotAlignedPositions(self):
        if not self.options.drawGridForNotAlignedPositions:
            return

        alignedPositions: List[XmapAlignedPair] = \
            [position for segment in self.alignment.segments for position in segment.alignedPositions] \
            + [pair for pair in self.benchmarkAlignment.alignedPairs]

        x = [p for p in self.reference.positions
             if self.__isReferencePositionInScope(p) and self.__isNotAlignedReference(p, alignedPositions)]
        y = [p for p in self.query.positions
             if self.__isQueryPositionInScope(p) and self.__isNotAlignedQuery(p, alignedPositions)]
        self.__drawGrid(x, y, self.yMaxPlot, self.xMaxPlot, lineStyle=(0, (5, 15)), label="not aligned positions")

    @staticmethod
    def __isNotAlignedReference(position: int, alignedPositions: List[XmapAlignedPair]):
        return position not in [ap.reference.position for ap in alignedPositions]

    @staticmethod
    def __isNotAlignedQuery(position: int, alignedPositions: List[XmapAlignedPair]):
        return position not in [ap.query.position for ap in alignedPositions]

    def __isReferencePositionInScope(self, position: int):
        return self.xMinPlot <= position <= self.xMaxAxis

    def __isQueryPositionInScope(self, position: int):
        return self.yMinAxis <= position <= self.yMaxPlot