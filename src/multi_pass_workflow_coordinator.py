import os
from typing import List, Dict, Tuple

from src.alignment.aligner import Aligner, ChainBuilder
from src.alignment.alignment_results import AlignmentResultRow, AlignmentResults
from src.args import Args
from src.correlation.optical_map import OpticalMap
from src.correlation.peaks_selector import PeaksSelector
from src.correlation.sequence_generator import SequenceGenerator
from src.extensions.dispatcher import Dispatcher
from src.parsers.xmap_reader import XmapReader
from src.workflow_coordinator import _WorkflowCoordinator


class _MultiPassWorkflowCoordinator(_WorkflowCoordinator):
    def __init__(self,
                 args: Args,
                 primaryGenerator: SequenceGenerator,
                 secondaryGenerator: SequenceGenerator,
                 aligner: Aligner | ChainBuilder,
                 dispatcher: Dispatcher,
                 peaksSelector: PeaksSelector,
                 xmapReader: XmapReader):
        super().__init__(args, primaryGenerator, secondaryGenerator, aligner, dispatcher, peaksSelector)
        self.xmapReader = xmapReader

    def execute(self, referenceMaps: List[OpticalMap], queryMaps: List[OpticalMap]) -> List[AlignmentResultRow]:
        pass_rows: List[List[AlignmentResultRow]] = []
        inactive_queries: set[int] = set()

        first_pass = super().execute(referenceMaps, queryMaps)
        query_len_by_id = {q.moleculeId: q.length for q in queryMaps}
        pass_rows.append(first_pass)

        all_queries = {q.moleculeId for q in queryMaps}
        inactive_queries |= all_queries - {r.queryId for r in first_pass}

        if isinstance(self.aligner, ChainBuilder):
            pass_no = 2
            while True:
                latest_rows = pass_rows[-1] if pass_rows else []
                next_fragments: List[OpticalMap] = []

                for r in latest_rows:
                    if r.queryId in inactive_queries:
                        continue

                    rests = (getattr(r, 'rests', []) or [])

                    total_rest_bp = 0.0
                    for frag in rests:
                        pos = getattr(frag, 'positions', None)
                        if pos and len(pos) >= 2:
                            total_rest_bp += (max(pos) - min(pos))

                    for frag in rests:
                        pos = getattr(frag, 'positions', None)
                        if pos is not None and len(pos) >= 7:
                            next_fragments.append(frag)

                if not next_fragments:
                    break

                new_rows = super().execute(referenceMaps, next_fragments, isFirstPass=False)

                attempted = {f.moleculeId for f in next_fragments}
                successful = {r.queryId for r in new_rows}
                inactive_queries |= attempted - successful

                if not new_rows:
                    break

                new_rows = [r.setAlignedRest(True) for r in new_rows]
                pass_rows.append(new_rows)
                pass_no += 1

        else:
            pass_no = 2
            while True:
                next_fragments = self._get_next_pass_fragments(pass_rows, queryMaps, inactive_queries)
                if not next_fragments:
                    break

                new_rows = super().execute(referenceMaps, next_fragments)

                attempted = {f.moleculeId for f in next_fragments}
                successful = {r.queryId for r in new_rows}
                inactive_queries |= attempted - successful

                if not new_rows:
                    break

                new_rows = [r.setAlignedRest(True) for r in new_rows]
                pass_rows.append(new_rows)
                pass_no += 1

        if self.args.outputMode == 'split-alignments':
            all_rows_flat = [row for rows in pass_rows for row in rows]
            return all_rows_flat

        filtered_per_pass = [
            AlignmentResults.filterOutSubsequentAlignmentsForSingleQuery(rows)
            for rows in pass_rows
        ]

        if self.args.outputMode == 'separate':
            for i, rows in enumerate(filtered_per_pass[1:], start=1):
                self.saveAdditionalOutput(rows, i)
            return filtered_per_pass[0]

        flattened_filtered = [row for rows in filtered_per_pass for row in rows]
        joinedRows, separateRows = AlignmentResults.resolve(flattened_filtered, self.args.maxDifference)

        if self.args.outputMode == "best":
            joined_ids = {row.queryId for row in joinedRows}
            bestRows = [row for row in filtered_per_pass[0] if row.queryId not in joined_ids]
            return sorted(joinedRows + bestRows, key=lambda r: r.queryId)

        if self.args.outputMode == 'joined':
            self.saveAdditionalOutput(separateRows, 1)
            return joinedRows

        if self.args.outputMode == 'all':
            for i, rows in enumerate(filtered_per_pass, start=1):
                self.saveAdditionalOutput(rows, i)
            return joinedRows

    def _get_next_pass_fragments(
            self,
            pass_rows: List[List[AlignmentResultRow]],
            queryMaps: List[OpticalMap],
            inactive_queries: set[int]) -> List[OpticalMap]:

        all_rows_flat = [row for rows in pass_rows for row in rows]
        rows_by_query: Dict[int, List[AlignmentResultRow]] = {}
        for r in all_rows_flat:
            rows_by_query.setdefault(r.queryId, []).append(r)

        query_by_id = {q.moleculeId: q for q in queryMaps}
        result_fragments: List[OpticalMap] = []

        for qid, rows in rows_by_query.items():
            if qid in inactive_queries:
                continue

            query = query_by_id.get(qid)
            if not query:
                continue

            if self._coverage_over_threshold(rows, query):
                continue

            fragments = self._intersection_of_unaligned(rows, query)
            if not fragments:
                continue

            for (start_idx, end_idx) in fragments:
                positions = query.positions[start_idx:end_idx + 1]
                if len(positions) >= 7:
                    result_fragments.append(
                        OpticalMap(qid, query.length, positions, shift=start_idx)
                    )
        return result_fragments


    def _coverage_over_threshold(self, rows: List[AlignmentResultRow], query: OpticalMap) -> bool:
        intervals = sorted(
            (min(r.queryStartPosition, r.queryEndPosition), max(r.queryStartPosition, r.queryEndPosition))
            for r in rows
        )
        merged: List[Tuple[float, float]] = []
        for s, e in intervals:
            if not merged or s > merged[-1][1]:
                merged.append([s, e])
            else:
                merged[-1][1] = max(merged[-1][1], e)
        coverage = sum(e - s for s, e in merged)
        return coverage > 0.8 * query.length
    
    def _intersection_of_unaligned(
            self,
            rows: List[AlignmentResultRow],
            query: OpticalMap) -> List[Tuple[int, int]]:

        complements: List[List[Tuple[int, int]]] = []
        for r in rows:
            frags = r.getUnalignedFragments([query])
            print(frags)
            row_intervals = []
            for f in frags:
                start_idx = f.shift
                end_idx = f.shift + len(f.positions) - 1
                row_intervals.append((start_idx, end_idx))
            if not row_intervals:
                return []
            complements.append(row_intervals)

        intersection = complements[0]
        for other in complements[1:]:
            new_intersection: List[Tuple[int, int]] = []
            for s1, e1 in intersection:
                for s2, e2 in other:
                    s = max(s1, s2)
                    e = min(e1, e2)
                    if s <= e:
                        new_intersection.append((s, e))
            if not new_intersection:
                return []

            new_intersection.sort()
            merged: List[Tuple[int, int]] = []
            for s, e in new_intersection:
                if not merged or s > merged[-1][1] + 1:
                    merged.append([s, e])
                else:
                    merged[-1][1] = max(merged[-1][1], e)
            intersection = [tuple(m) for m in merged]

        return intersection

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

