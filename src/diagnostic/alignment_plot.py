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

from src.alignment.alignment_position import AlignedPair
from src.alignment.alignment_results import AlignmentResultRow
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import OpticalMap, PositionWithSiteId, InitialAlignment
from src.correlation.peak import Peak
from src.diagnostic.benchmark_alignment import BenchmarkAlignment, BenchmarkAlignedPair


@dataclass(frozen=True)
class Options:
    limitQueryToAlignedArea: bool = False
    drawGridForNotAlignedPositions: bool = True
    drawRemovedAlignedPositions: bool = True


class AlignmentPlot:
    def __init__(self, reference: OpticalMap, query: OpticalMap, alignment: AlignmentResultRow | BenchmarkAlignment,
                 correlation: InitialAlignment = None, benchmarkAlignment: BenchmarkAlignment = None, options: Options = None):
        self.reference = reference
        self.query = query
        self.alignment = alignment
        self.correlation = correlation
        self.benchmarkAlignment = benchmarkAlignment
        self.options = options or Options()
        self.create()

    def create(self):
        self._createFigure()
        self._setDimensions()
        self._plotReference()
        self._plotQuery()
        self._drawPrimaryPeak()
        self._plotBenchmarkAlignment()
        self._plotSegments()
        self._drawGridForNotAlignedPositions()
        self._drawRemovedAlignedPositions()
        self.axes.legend()

    def _createFigure(self):
        self.figure: Figure = pyplot.figure()
        self.axes: Axes = self.figure.add_axes((0, 0, 1, 1))
        self.axes.set_aspect("equal")
        self.axes.ticklabel_format(style='plain')
        formatter = EngFormatter(unit="bp")
        self.axes.xaxis.set_major_formatter(formatter)
        self.axes.yaxis.set_major_formatter(formatter)

    def _setDimensions(self):
        margin = self.alignment.queryLength / 10
        drawPeakPositions = \
            list(map(lambda s: s.peak.leftProminenceBasePosition, self.alignment.segments)) \
            + [self.correlation.maxPeak.leftProminenceBasePosition] \
                if self.correlation and self.correlation.maxPeak and hasattr(self.alignment, "segments") else []

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

    def _plotReference(self):
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
            self.axes.annotate(str(r.siteId), (r.position, self.yMinAxis),
                               textcoords="offset points",
                               xytext=(0, -10),
                               ha="center",
                               va="top",
                               rotation=90)

    def _plotQuery(self):
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
            self.axes.annotate(str(q.siteId), (self.xMinAxis, q.position),
                               textcoords="offset points",
                               xytext=(-15, 0),
                               ha="right",
                               va="center")

    def _plotSegments(self):
        groupedSegments = groupby(self.alignment.segments, lambda s: s.peak)
        colors = self.__getContrastingColors(len(self.alignment.segments))

        for (peakNumber, (peak, segments)), color in zip(enumerate(groupedSegments), colors):
            self.__drawPeak(color, peak)
            for segmentNumber, segment in enumerate(segments):
                x = list(map(lambda p: p.reference.position, segment.alignedPositions))
                y = list(map(lambda p: self.__absoluteQueryPosition(p), segment.alignedPositions))
                self.__plotSegment(color, peak, peakNumber, segment, segmentNumber, x, y)
                self.__drawGrid(x, y)

    def __drawPeak(self, color, peak: Peak):
        peakRectangle = Rectangle((peak.leftProminenceBasePosition, self.yMinBorder - 9999999999.),
                                  peak.width,
                                  self.yMaxAxis + 99999999999999.,
                                  color=color,
                                  alpha=0.2,
                                  angle=self.__drawPeakAngle,
                                  rotation_point=self.__drawPeakRotationPoint(peak))
        self.axes.add_patch(peakRectangle)
        peakRectangle.set_clip_path(self.plotAreaMask)

    def _drawPrimaryPeak(self):
        if not self.correlation and self.correlation.maxPeak:
            return

        peak = self.correlation.maxPeak
        x = [peak.position, self.xMaxPlot]
        y = [0, self.xMaxPlot - peak.position]

        def reverseY():
            return [self.query.length, self.query.length - (self.xMaxPlot - peak.position)]

        self.axes.plot(x, reverseY() if self.alignment.reverseStrand else y,
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
                                  angle=self.__drawPeakAngle,
                                  rotation_point=self.__drawPeakRotationPoint(peak))
        self.axes.add_patch(peakRectangle)
        peakRectangle.set_clip_path(self.plotAreaMask)

    def _plotBenchmarkAlignment(self):
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
                       label=f" peak {peakNumber + 1} (height: {peak.height:.2f}), segment {segmentNumber + 1} "
                             f"({len(segment.alignedPositions)} pairs, score: {segment.segmentScore:.1f})",
                       marker="+",
                       markersize=16,
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

    def _drawGridForNotAlignedPositions(self):
        if not self.options.drawGridForNotAlignedPositions or not hasattr(self.alignment, "segments"):
            return

        alignedPositions: List[BenchmarkAlignedPair] = \
            [position for segment in self.alignment.segments for position in segment.alignedPositions] \
            + [pair for pair in self.benchmarkAlignment.alignedPairs]

        x = [p for p in self.reference.positions if
             self.__isReferencePositionInScope(p) and self.__isNotAlignedReference(p, alignedPositions)]
        y = [p for p in self.query.positions if
             self.__isQueryPositionInScope(p) and self.__isNotAlignedQuery(p, alignedPositions)]
        self.__drawGrid(x, y, self.yMaxPlot, self.xMaxPlot, lineStyle=(0, (15, 3)), label="not aligned positions")

    def _drawRemovedAlignedPositions(self):
        if not self.options.drawRemovedAlignedPositions or not hasattr(self.alignment, "segments"):
            return

        segmentsAlignedPositions: List[BenchmarkAlignedPair] = \
            [position for segment in self.alignment.segments for position in segment.alignedPositions]

        allAlignedPositions: List[AlignedPair] = \
            list(set(position for segment in self.alignment.segments for position in segment.allPeakPositions if
                     isinstance(position, AlignedPair)))

        removedAlignedPositions = sorted([p for p in allAlignedPositions if p not in segmentsAlignedPositions])

        self.axes.scatter(
            [p.reference.position for p in removedAlignedPositions],
            [self.__absoluteQueryPosition(p) for p in removedAlignedPositions],
            label=f"aligned pairs removed from segments ({len(removedAlignedPositions)} pairs)",
            marker="|",
            s=16 ** 2,
            c="orange")

        self.__annotateSources(removedAlignedPositions)

    def __annotateSources(self, positions):
        for position, groups in groupby(positions):
            samePositionsFromDifferentSources = list(groups)
            position = samePositionsFromDifferentSources[0]
            sources = map(str, sorted(map(lambda p: p.source, samePositionsFromDifferentSources)))
            self.axes.annotate(f"peak:{','.join(sources)}",
                               (position.reference.position, self.__absoluteQueryPosition(position)))

    @staticmethod
    def __isNotAlignedReference(position: int, alignedPositions: List[BenchmarkAlignedPair]):
        return position not in [ap.reference.position for ap in alignedPositions]

    @staticmethod
    def __isNotAlignedQuery(position: int, alignedPositions: List[BenchmarkAlignedPair]):
        return position not in [ap.query.position for ap in alignedPositions]

    def __isReferencePositionInScope(self, position: int):
        return self.xMinPlot <= position <= self.xMaxAxis

    def __isQueryPositionInScope(self, position: int):
        return self.yMinAxis <= position <= self.yMaxPlot

    def __absoluteQueryPosition(self, p: BenchmarkAlignedPair | AlignedPair):
        return self.alignment.queryLength - p.query.position if self.alignment.reverseStrand else p.query.position

    @property
    def __drawPeakAngle(self):
        return 45 if self.alignment.reverseStrand else -45.

    def __drawPeakRotationPoint(self, peak: Peak):
        return peak.position, self.query.length if self.alignment.reverseStrand else 0


class BenchmarkAlignmentPlot(AlignmentPlot):
    def __init__(self, reference: OpticalMap, query: OpticalMap, alignment: BenchmarkAlignment, options: Options = None):
        super().__init__(reference, query, alignment, None, None, options)

    def create(self):
        self._createFigure()
        self._setDimensions()
        self._plotReference()
        self._plotQuery()
        self._plotAlignment()
        self.axes.legend()

    def _plotAlignment(self):
        x, y = list(zip(*map(lambda p: (p.reference.position, p.query.position), self.alignment.alignedPairs)))
        self.axes.plot(x, y,
                       label=f"({len(self.alignment.alignedPairs)} pairs, "
                             f"confidence: {self.alignment.confidence})",
                       color="red",
                       markeredgecolor="red",
                       fillstyle="none",
                       marker="o",
                       markersize=8,
                       linewidth=2)
