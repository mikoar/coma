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
        self.options = options or Options()
        self.__createFigure()
        self.__setDimensions(alignment, benchmarkAlignment, query, correlation)
        self.__plotReference(reference)
        self.__plotQuery(query)
        self.__drawPrimaryPeak(correlation)
        self.__plotBenchmarkAlignment(benchmarkAlignment)
        self.__plotSegments(alignment)
        self.__drawGridForNotAlignedPositions(reference, query, alignment, benchmarkAlignment)
        self.axes.legend()

    def __createFigure(self):
        self.figure: Figure = pyplot.figure()
        self.axes: Axes = self.figure.add_axes([0, 0, 1, 1])
        self.axes.set_aspect("equal")
        self.axes.ticklabel_format(style='plain')
        formatter = EngFormatter(unit="bp")
        self.axes.xaxis.set_major_formatter(formatter)
        self.axes.yaxis.set_major_formatter(formatter)

    def __setDimensions(self, alignment: AlignmentResultRow, benchmarkAlignment: XmapAlignment, query: OpticalMap,
                        correlation: InitialAlignment):
        offset = alignment.queryLength / 10
        drawPeakPositionsWithOffset = \
            list(map(lambda s: s.peak.leftProminenceBasePosition - offset, alignment.segments)) \
            + [correlation.maxPeak.leftProminenceBasePosition - offset]

        self.overlapsWithBenchmark = (alignment.referenceStartPosition <= benchmarkAlignment.referenceStartPosition
                                      <= alignment.referenceEndPosition or alignment.referenceStartPosition
                                      <= benchmarkAlignment.referenceEndPosition <= alignment.referenceEndPosition) \
            if benchmarkAlignment else False

        alignmentStartPositionsWithOffset = [alignment.referenceStartPosition - offset]
        alignmentReferenceEndPositions = [alignment.referenceEndPosition]

        self.yMin = (alignment.queryStartPosition if self.options.limitQueryToAlignedArea else 0) - offset
        self.yMax = alignment.queryEndPosition if self.options.limitQueryToAlignedArea else query.length

        self.yMinWithOffset = self.yMin - offset
        self.yMaxWithOffset = self.yMax + offset

        if self.overlapsWithBenchmark:
            alignmentStartPositionsWithOffset.append(benchmarkAlignment.referenceStartPosition - offset)
            alignmentReferenceEndPositions.append(benchmarkAlignment.referenceEndPosition)

        self.xMin = min(alignmentStartPositionsWithOffset + drawPeakPositionsWithOffset)
        self.xMax = max(alignmentReferenceEndPositions)

        self.xMinWithOffset = self.xMin - offset
        self.xMaxWithOffset = self.xMax + offset
        self.axes.set_xlim(self.xMinWithOffset, self.xMaxWithOffset)
        self.axes.set_ylim(self.yMinWithOffset, self.yMaxWithOffset)
        self.plotAreaMask = Rectangle((self.xMin + offset, self.yMin + offset),
                                      self.xMax - self.xMin - offset,
                                      self.yMax - self.yMin - offset,
                                      facecolor='none', edgecolor='none')
        self.axes.add_patch(self.plotAreaMask)
        self.__setFigureSize()

    def __setFigureSize(self):
        scale = 45_000
        minSize = 10.
        maxSize = 60.
        xSize = (self.xMaxWithOffset - self.xMinWithOffset) / scale
        ySize = (self.yMaxWithOffset - self.yMinWithOffset) / scale
        xSizeClamped = min(max([xSize, minSize]), maxSize)
        ySizeClamped = min(max([ySize, minSize]), maxSize)
        self.figure.set_size_inches((xSizeClamped, ySizeClamped))

    def __plotReference(self, reference: OpticalMap):
        self.axes.set_xlabel(f"Reference {reference.moleculeId}")
        refLabelsInScope = [p for p in reference.getPositionsWithSiteIds()
                            if self.__isReferencePositionInScope(p.position)]
        self.axes.plot([r.position for r in refLabelsInScope],
                       np.repeat(self.yMin, len(refLabelsInScope)),
                       marker="|",
                       markersize=16,
                       markeredgewidth="2", linewidth=16, markeredgecolor="black",
                       color="yellow")

        for r in self.__skipDensePositions(refLabelsInScope):
            self.axes.annotate(r.siteId, (r.position, self.yMin),
                               textcoords="offset points",
                               xytext=(0, -10),
                               ha="center",
                               va="top",
                               rotation=90)

    def __plotQuery(self, query: OpticalMap):
        self.axes.set_ylabel(f"Query {query.moleculeId}")
        queryLabelsInScope = [p for p in query.getPositionsWithSiteIds()
                              if not self.options.limitQueryToAlignedArea or self.__isQueryPositionInScope(p.position)]

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
            self.__drawPeak(color, peak)
            for segmentNumber, segment in enumerate(segments):
                x = list(map(lambda p: p.reference.position, segment.alignedPositions))
                y = list(map(lambda p: p.query.position, segment.alignedPositions))
                self.__plotSegment(color, peak, peakNumber, segment, segmentNumber, x, y)
                self.__drawGrid(x, y)

    def __drawPeak(self, color, peak: Peak):
        peakRectangle = Rectangle((peak.leftProminenceBasePosition, self.yMinWithOffset - 9999999999.),
                                  peak.width,
                                  self.yMaxWithOffset + 99999999999999.,
                                  color=color,
                                  alpha=0.2,
                                  angle=-45.,
                                  rotation_point=(peak.position, 0))  # type:ignore
        self.axes.add_patch(peakRectangle)
        peakRectangle.set_clip_path(self.plotAreaMask)

    def __drawPrimaryPeak(self, correlation: InitialAlignment):
        peak = correlation.maxPeak
        self.axes.plot([peak.position, self.xMax], [0, self.xMax - peak.position],
                       linestyle="dashdot",
                       marker=None,
                       color="black")

        peakRectangle = Rectangle((peak.leftProminenceBasePosition, self.yMinWithOffset - 9999999999.),
                                  peak.width,
                                  self.yMaxWithOffset + 99999999999999.,
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

    def __plotBenchmarkAlignment(self, benchmarkAlignment: XmapAlignment):
        if not benchmarkAlignment:
            return

        x, y = list(zip(*map(lambda p: (p.reference.position, p.query.position), benchmarkAlignment.alignedPairs)))
        self.axes.plot(x, y,
                       label=f"benchmark ({len(benchmarkAlignment.alignedPairs)} pairs, "
                             f"confidence: {benchmarkAlignment.confidence})",
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

    def __drawGrid(self, x, y, xHeights=None, yHeights=None, lineStyle: str | Tuple = "--"):
        self.axes.vlines(x, self.yMin, xHeights or y, linestyles=lineStyle, colors="gray", linewidth=0.5)
        self.axes.hlines(y, self.xMin, yHeights or x, linestyles=lineStyle, colors="gray", linewidth=0.5)

    def __skipDensePositions(self, positions: List[PositionWithSiteId]):
        minDistance = (self.xMax - self.xMin) / 150
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

    def __drawGridForNotAlignedPositions(self, reference: OpticalMap, query: OpticalMap, alignment: AlignmentResultRow,
                                         benchmarkAlignment: XmapAlignment):
        if not self.options.drawGridForNotAlignedPositions:
            return

        alignedPositions: List[XmapAlignedPair] = \
            [position for segment in alignment.segments for position in segment.alignedPositions] \
            + [pair for pair in benchmarkAlignment.alignedPairs]

        x = [p for p in reference.positions
             if self.__isReferencePositionInScope(p) and self.__isNotAlignedReference(p, alignedPositions)]
        y = [p for p in query.positions
             if self.__isQueryPositionInScope(p) and self.__isNotAlignedQuery(p, alignedPositions)]
        self.__drawGrid(x, y, self.yMax, self.xMax, lineStyle=(0, (5, 15)))

    @staticmethod
    def __isNotAlignedReference(position: int, alignedPositions: List[XmapAlignedPair]):
        return position not in [ap.reference.position for ap in alignedPositions]

    @staticmethod
    def __isNotAlignedQuery(position: int, alignedPositions: List[XmapAlignedPair]):
        return position not in [ap.query.position for ap in alignedPositions]

    def __isReferencePositionInScope(self, position: int):
        return self.xMin <= position <= self.xMaxWithOffset

    def __isQueryPositionInScope(self, position: int):
        return self.yMin <= position <= self.yMax
