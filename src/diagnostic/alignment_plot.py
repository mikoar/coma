from itertools import groupby

import numpy as np
from matplotlib import pyplot
from matplotlib.axes import Axes
from matplotlib.figure import Figure

from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.optical_map import OpticalMap


class AlignmentPlot:
    def __init__(self, reference: OpticalMap, query: OpticalMap, alignment: AlignmentResultRow):
        self.figure: Figure = pyplot.figure(figsize=(24, 16))
        self.axes: Axes = self.figure.add_axes([0, 0, 1, 1])
        self.__setLimits(alignment)
        self.__plotReference(alignment, reference)
        self.__plotQuery(alignment, query)
        self.__plotSegments(alignment)
        self.axes.legend()

    def __setLimits(self, alignment: AlignmentResultRow):
        offset = (alignment.queryEndPosition - alignment.queryStartPosition) / 10
        self.yMin = - offset
        peakPositions = list(map(lambda s: s.peak, alignment.segments))
        self.xMin = min([alignment.referenceStartPosition - offset] + peakPositions)
        self.axes.set_xlim(self.xMin - offset, alignment.referenceEndPosition + offset)
        self.axes.set_ylim(self.yMin - offset, alignment.queryLength + offset)
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
                       color=self.__getColors(2)[0])

        for r in referenceLabelsInScope:
            self.axes.annotate(r.siteId, (r.position, self.yMin),
                               textcoords="offset points",
                               xytext=(0, -20),
                               ha="center",
                               rotation=90)

    def __plotQuery(self, alignment: AlignmentResultRow, query: OpticalMap):
        self.axes.set_ylabel(f"Query {query.moleculeId}")
        queryLabelsInScope = [p for p in query.getPositionsWithSiteIds() if
                              self.yMin <= p.position <= alignment.queryEndPosition]
        self.axes.plot(np.repeat(self.xMin, len(queryLabelsInScope)), [q.position for q in queryLabelsInScope],
                       marker="_",
                       markersize=16,
                       markeredgewidth="2", linewidth=16, markeredgecolor="black",
                       color=self.__getColors(2)[1])

        for q in queryLabelsInScope:
            self.axes.annotate(q.siteId, (self.xMin, q.position), textcoords="offset points", xytext=(-20, 0),
                               va="center")

    def __plotSegments(self, alignment: AlignmentResultRow):
        groupedSegments = groupby(alignment.segments, lambda s: s.peak)
        colors = self.__getColors(len(alignment.segments))

        for (peakNumber, (_, segments)), color in zip(enumerate(groupedSegments), colors):
            for segmentNumber, segment in enumerate(segments):
                x = list(map(lambda p: p.reference.position, segment.alignedPositions))
                y = list(map(lambda p: p.query.position, segment.alignedPositions))
                self.axes.plot(x, y,
                               label=f" peak {peakNumber + 1} segment {segmentNumber + 1}",
                               marker="+",
                               markersize=8,
                               markeredgecolor="black",
                               linewidth=2,
                               color=color)

                self.axes.vlines(x, self.yMin, y, linestyles="--", colors="black", alpha=0.5, linewidth=1)
                self.axes.hlines(y, self.xMin, x, linestyles="--", colors="black", alpha=0.5, linewidth=1)

    @staticmethod
    def __getColors(count):
        colorMap = pyplot.get_cmap("plasma")
        return colorMap(np.linspace(0, 1, count))
