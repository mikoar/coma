import os
from typing import List

from alignment.aligner import Aligner
from alignment.alignment_results import AlignmentResultRow, AlignmentResults
from args import Args
from correlation.optical_map import OpticalMap
from correlation.peaks_selector import PeaksSelector
from correlation.sequence_generator import SequenceGenerator
from extensions.dispatcher import Dispatcher
from parsers.xmap_reader import XmapReader
from workflow_coordinator import _WorkflowCoordinator


class _MultiPassWorkflowCoordinator(_WorkflowCoordinator):
    def __init__(self,
                 args: Args,
                 primaryGenerator: SequenceGenerator,
                 secondaryGenerator: SequenceGenerator,
                 aligner: Aligner,
                 dispatcher: Dispatcher,
                 peaksSelector: PeaksSelector,
                 xmapReader: XmapReader):
        super().__init__(args, primaryGenerator, secondaryGenerator, aligner, dispatcher, peaksSelector)
        self.xmapReader = xmapReader

    def execute(self, referenceMaps: List[OpticalMap], queryMaps: List[OpticalMap]) -> List[AlignmentResultRow]:
        alignmentResultRows = super().execute(referenceMaps, queryMaps)
        alignmentResultRowsSecondPass = self.getSecondPassAlignmentRows(alignmentResultRows, queryMaps, referenceMaps)

        if self.args.outputMode == "best":
            alignmentResultRows = alignmentResultRows + alignmentResultRowsSecondPass

        filteredFirstPassRows = AlignmentResults.filterOutSubsequentAlignmentsForSingleQuery(alignmentResultRows)
        filteredSecondPassRows = \
            AlignmentResults.filterOutSubsequentAlignmentsForSingleQuery(alignmentResultRowsSecondPass)

        if self.args.outputMode == 'separate':
            self.saveAdditionalOutput(filteredSecondPassRows)
            return filteredFirstPassRows

        joinedRows, separateRows = AlignmentResults.resolve(
            filteredFirstPassRows + filteredSecondPassRows,
            self.args.maxDifference)

        if self.args.outputMode == 'best':
            joinedIds = [row.queryId for row in joinedRows]
            bestRows = [row for row in filteredFirstPassRows if row.queryId not in joinedIds]
            bestAndJoinedRows = sorted(joinedRows + bestRows, key=lambda r: r.queryId)
            return bestAndJoinedRows

        if self.args.outputMode == 'joined':
            self.saveAdditionalOutput(separateRows, 1)
            return joinedRows

        if self.args.outputMode == 'all':
            self.saveAdditionalOutput(filteredFirstPassRows, 1)
            self.saveAdditionalOutput(filteredSecondPassRows, 2)
            return joinedRows

    def getSecondPassAlignmentRows(self, alignmentResultRows, queryMaps, referenceMaps):
        unalignedFragmentsLists = \
            [alignmentResultRow.getUnalignedFragments(queryMaps) for alignmentResultRow in alignmentResultRows]
        unalignedFragments = [item for row in unalignedFragmentsLists for item in row]
        alignmentResultRowsSecondPass = super().execute(referenceMaps, unalignedFragments)
        alignmentResultRowsSecondPass = [alignmentResultRowRest.setAlignedRest(True) for alignmentResultRowRest in
                                         alignmentResultRowsSecondPass]
        return alignmentResultRowsSecondPass

    def saveAdditionalOutput(
            self,
            rowsWithoutSubsequentAlignmentsForSingleQueryRest: List[AlignmentResultRow],
            fileNumber: int):
        restResult = AlignmentResults(
            self.args.referenceFile.name,
            self.args.queryFile.name,
            rowsWithoutSubsequentAlignmentsForSingleQueryRest)

        self.xmapReader.writeAlignments(self.createAdditionalOutputFile(fileNumber), restResult, self.args)

    def createAdditionalOutputFile(self, number: int):
        return open("{0}_{2}{1}".format(*os.path.splitext(self.args.outputFile.name) + (number,)), mode='w',
                    encoding=self.args.outputFile.encoding)
