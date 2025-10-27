from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from enum import Enum
from typing import Optional, Dict, Any, ClassVar, List, Tuple
from collections import defaultdict
import math
import re

from src.alignment.alignment_position import AlignedPair, NotAlignedPosition
from src.alignment.segment_with_resolved_conflicts import AlignmentSegmentsWithResolvedConflicts
from src.alignment.segments import AlignmentSegment
from src.correlation.optical_map import OpticalMap
from src.diagnostic.benchmark_alignment import BenchmarkAlignment


class HitEnum(Enum):
    MATCH = "M"
    DELETION = "D"
    INSERTION = "I"

@dataclass
class SmapResultRow:

    SV_PROXIMITY_THR: int = 50_000

    COLUMN_NAMES = [
        "SmapEntryID","QryContigID","RefcontigID1","RefcontigID2",
        "QryStartPos","QryEndPos","RefStartPos","RefEndPos",
        "Confidence","Type","XmapID1","XmapID2","LinkID",
        "QryStartIdx","QryEndIdx","RefStartIdx","RefEndIdx",
        "RawConfidence","RawConfidenceLeft","RawConfidenceRight","RawConfidenceCenter",
        "SVsize","SVfreq","Orientation","VAF",
    ]
    COLUMN_TYPES = [
        "int","int","int","int",
        "float","float","float","float",
        "float","string","int","int","int",
        "int","int","int","int",
        "float","float","float","float",
        "float","float","string","float",
    ]

    def __init__(self) -> None:
        self.data: Dict[str, Any] = {k: None for k in self.COLUMN_NAMES}
        self.data["RawConfidence"] = -1.0
        self.data["RawConfidenceCenter"] = -1.0
        self.data["RawConfidenceLeft"] = 50.0
        self.data["RawConfidenceRight"] = 50.0
        self.data["SVfreq"] = -1.0
        self.data["VAF"] = -1.0

    def to_smap_dict(self) -> Dict[str, Any]:
        return dict(self.data)

    @classmethod
    def from_alignment_group(cls, qid: int, group_rows: List[Any]) -> List["SmapResultRow"]:
        n = len(group_rows)
        if n < 2:
            return []

        out: List[SmapResultRow] = []
        used: set[int] = set()

        for i in range(n - 3):
            if {i, i+1, i+2, i+3} & used:
                continue
            L, A, B, R = group_rows[i:i+4]
            if cls.is_full_duplication_window(L, A, B, R):
                row = cls.build_entry_pair(qid, A, B, sv_type="duplication", orientation_str=cls.ori_char(A))
                row.data["LinkID"] = -1
                out.append(row)
                used.update({i, i+1, i+2, i+3})

        inv_rows_idx: List[int] = []
        for i in range(n - 1):
            if {i, i+1} & used:
                continue
            A, B = group_rows[i], group_rows[i+1]
            if cls.is_inversion_pair(A, B):
                svt = "Inversion" if len(inv_rows_idx) == 0 else "Inversion_paired"
                row = cls.build_entry_pair(qid, A, B, sv_type=svt, orientation_str=None)
                row.data["LinkID"] = -1
                out.append(row)
                inv_rows_idx.append(len(out) - 1)

        if len(inv_rows_idx) >= 2:
            i1, i2 = inv_rows_idx[0], inv_rows_idx[1]
            setattr(out[i1], "_link_peer_idx", i2)
            setattr(out[i2], "_link_peer_idx", i1)

        for i in range(n - 1):
            if {i, i+1} & used:
                continue
            A, B = group_rows[i], group_rows[i+1]
            if cls.is_duplication_split_pair(A, B):
                row = cls.build_entry_pair(qid, A, B, sv_type="duplication_split", orientation_str=cls.ori_char(A))
                row.data["LinkID"] = -1
                out.append(row)

        for i in range(n - 1):
            if {i, i+1} & used:
                continue
            A, B = group_rows[i], group_rows[i+1]
            t = cls.translocation_type(A, B)
            if t is None:
                continue
            orient = cls.pair_orientation(A, B)
            row = cls.build_entry_pair(qid, A, B, sv_type=t, orientation_str=orient)
            row.data["LinkID"] = -1
            out.append(row)

        return out


    @staticmethod
    def ori_char(row: Any) -> str:
        if hasattr(row, "orientation") and row.orientation in ("+","-"):
            return row.orientation
        return "-" if getattr(row, "reverseStrand", False) else "+"

    @staticmethod
    def pair_orientation(a: Any, b: Any) -> str:
        return f"{SmapResultRow.ori_char(a)}/{SmapResultRow.ori_char(b)}"

    @staticmethod
    def rid(row: Any) -> Optional[int]:
        v = getattr(row, "referenceId", None)
        try: return int(v) if v is not None else None
        except Exception: return None

    @staticmethod
    def get_xmap_id(row: Any) -> Optional[int]:
        for cand in ("xmapId","XmapEntryID","xmap_entry_id","entryId","id"):
            v = getattr(row, cand, None)
            if v is not None:
                try: return int(v)
                except Exception: pass
        return None

    @classmethod
    def iv(cls, lo: float, hi: float) -> Tuple[float,float]:
        return (min(lo,hi), max(lo,hi))

    @classmethod
    def ref_iv(cls, r: Any) -> Tuple[float,float]:
        return cls.iv(float(r.referenceStartPosition), float(r.referenceEndPosition))

    @classmethod
    def qry_iv(cls, r: Any) -> Tuple[float,float]:
        return cls.iv(float(r.queryStartPosition), float(r.queryEndPosition))

    @classmethod
    def span_both_axes(cls, a: Any, b: Any) -> Tuple[float,float,float,float]:
        a_qlo,a_qhi = cls.qry_iv(a); b_qlo,b_qhi = cls.qry_iv(b)
        a_rlo,a_rhi = cls.ref_iv(a); b_rlo,b_rhi = cls.ref_iv(b)
        return (min(a_qlo,b_qlo), max(a_qhi,b_qhi),
                min(a_rlo,b_rlo), max(a_rhi,b_rhi))

    @classmethod
    def min_conf(cls, *rows: Any) -> Optional[float]:
        vals = []
        for r in rows:
            c = getattr(r, "confidence", None)
            try:
                c = float(c)
                if not math.isnan(c):
                    vals.append(c)
            except Exception:
                pass
        return float(min(vals)) if vals else None

    @classmethod
    def ends_proximity(cls, a: Any, b: Any) -> float:
        a_lo,a_hi = cls.ref_iv(a); b_lo,b_hi = cls.ref_iv(b)
        return min(abs(a_lo-b_lo), abs(a_lo-b_hi), abs(a_hi-b_lo), abs(a_hi-b_hi))

    @classmethod
    def same_ref(cls, a: Any, b: Any) -> bool:
        ra, rb = cls.rid(a), cls.rid(b)
        return (ra is not None) and (ra == rb)

    @classmethod
    def same_ref_same_ori(cls, a: Any, b: Any) -> bool:
        return cls.same_ref(a,b) and (cls.ori_char(a) == cls.ori_char(b))

    @classmethod
    def overlap_or_close_on_ref(cls, a: Any, b: Any) -> bool:
        a_lo,a_hi = cls.ref_iv(a); b_lo,b_hi = cls.ref_iv(b)
        if min(a_hi,b_hi) - max(a_lo,b_lo) >= 0:
            return True
        return min(abs(a_hi - b_lo), abs(b_hi - a_lo)) <= float(cls.SV_PROXIMITY_THR)

    @classmethod
    def far_on_ref(cls, a: Any, b: Any) -> bool:
        a_lo,a_hi = cls.ref_iv(a); b_lo,b_hi = cls.ref_iv(b)
        if min(a_hi,b_hi) - max(a_lo,b_lo) >= 0:
            return False
        gap = min(abs(a_hi - b_lo), abs(b_hi - a_lo))
        return gap > float(cls.SV_PROXIMITY_THR)

    @classmethod
    def adjacent_on_qry(cls, a: Any, b: Any) -> bool:
        a_lo,a_hi = cls.qry_iv(a); b_lo,_ = cls.qry_iv(b)
        gap = max(0.0, b_lo - a_hi)
        return gap <= float(cls.SV_PROXIMITY_THR)

    @staticmethod
    def parse_alignment_pairs(row: Any) -> list[tuple[int, int]]:
        s = getattr(row, "alignment", None) or getattr(row, "Alignment", None) \
            or getattr(row, "alignmentPairs", None) or getattr(row, "align_pairs", None)
        if not s or not isinstance(s, str):
            return []
        pairs = []
        for m in re.finditer(r"\(\s*(\d+)\s*,\s*(\d+)\s*\)", s):
            r_idx = int(m.group(1))
            q_idx = int(m.group(2))
            pairs.append((r_idx, q_idx))
        return pairs

    @classmethod
    def edge_label_indices(cls, A: Any, B: Any) -> tuple[Optional[int], Optional[int], Optional[int], Optional[int]]:
        a_pairs = cls.parse_alignment_pairs(A)
        b_pairs = cls.parse_alignment_pairs(B)
        if not a_pairs or not b_pairs:
            return (None, None, None, None)

        a_ref_end, a_qry_end = a_pairs[-1]
        b_ref_start, b_qry_start = b_pairs[0]

        QryStartIdx = a_qry_end
        QryEndIdx   = b_qry_start
        RefStartIdx = a_ref_end
        RefEndIdx   = b_ref_start
        return (QryStartIdx, QryEndIdx, RefStartIdx, RefEndIdx)

    @staticmethod
    def labels_minmax_from_pairs(row: Any) -> Optional[Tuple[int, int, int, int]]:
        pairs = getattr(row, "alignedPairs", None)
        if not pairs:
            return None
        try:
            ref_ids = [p.reference.siteId for p in pairs]
            qry_ids = [p.query.siteId for p in pairs]
            return (min(qry_ids), max(qry_ids), min(ref_ids), max(ref_ids))
        except Exception:
            return None

    _ALIGN_RE = re.compile(r"\((\d+),(\d+)\)")

    @classmethod
    def labels_minmax_from_alignment_str(cls, row: Any) -> Optional[Tuple[int, int, int, int]]:
        s = getattr(row, "Alignment", None)
        if not s:
            return None
        ref_ids: List[int] = []
        qry_ids: List[int] = []
        for m in cls._ALIGN_RE.finditer(s):
            ref_ids.append(int(m.group(1)))
            qry_ids.append(int(m.group(2)))
        if not ref_ids or not qry_ids:
            return None
        return (min(qry_ids), max(qry_ids), min(ref_ids), max(ref_ids))

    @classmethod
    def labels_minmax(cls, row: Any) -> Optional[Tuple[int, int, int, int]]:
        mm = cls.labels_minmax_from_pairs(row)
        if mm is not None:
            return mm
        return cls.labels_minmax_from_alignment_str(row)

    @classmethod
    def fill_indices_from_rows(cls, A: Any, B: Any, out_row: "SmapResultRow") -> None:
        mmA = cls.labels_minmax(A)
        mmB = cls.labels_minmax(B)
        if mmA is None and mmB is None:
            return

        qmins, qmaxs, rmins, rmaxs = [], [], [], []
        if mmA:
            qa_min, qa_max, ra_min, ra_max = mmA
            qmins.append(qa_min); qmaxs.append(qa_max)
            rmins.append(ra_min); rmaxs.append(ra_max)
        if mmB:
            qb_min, qb_max, rb_min, rb_max = mmB
            qmins.append(qb_min); qmaxs.append(qb_max)
            rmins.append(rb_min); rmaxs.append(rb_max)

        if qmins and qmaxs:
            out_row.data["QryStartIdx"] = int(min(qmins))
            out_row.data["QryEndIdx"]   = int(max(qmaxs))
        if rmins and rmaxs:
            out_row.data["RefStartIdx"] = int(min(rmins))
            out_row.data["RefEndIdx"]   = int(max(rmaxs))

    @classmethod
    def is_duplication_split_pair(cls, A: Any, B: Any) -> bool:
        return cls.same_ref_same_ori(A,B) and cls.overlap_or_close_on_ref(A,B) and cls.adjacent_on_qry(A,B)

    @classmethod
    def is_full_duplication_window(cls, L: Any, A: Any, B: Any, R: Any) -> bool:
        if not cls.is_duplication_split_pair(A,B):
            return False
        rids = {cls.rid(x) for x in (L,A,B,R)}
        if len(rids) != 1 or (None in rids):
            return False
        if not (cls.ori_char(L) == cls.ori_char(A) == cls.ori_char(B) == cls.ori_char(R)):
            return False
        if not (cls.adjacent_on_qry(L,A) and cls.adjacent_on_qry(A,B) and cls.adjacent_on_qry(B,R)):
            return False
        if cls.ends_proximity(L,A) > float(cls.SV_PROXIMITY_THR): return False
        if cls.ends_proximity(B,R) > float(cls.SV_PROXIMITY_THR): return False
        return True

    @classmethod
    def is_inversion_pair(cls, A: Any, B: Any) -> bool:
        return cls.same_ref(A,B) and (cls.ori_char(A) != cls.ori_char(B)) \
               and (cls.ends_proximity(A,B) <= float(cls.SV_PROXIMITY_THR)) and cls.adjacent_on_qry(A,B)

    @classmethod
    def translocation_type(cls, A: Any, B: Any) -> Optional[str]:
        if not cls.adjacent_on_qry(A,B):
            return None
        if not cls.same_ref(A,B):
            return "translocation_interchr"
        if cls.far_on_ref(A,B):
            return "translocation_intrachr"
        return None

    @classmethod
    def build_entry_pair(cls,
                          qid: int,
                          A: Any,
                          B: Any,
                          sv_type: str,
                          orientation_str: Optional[str]) -> "SmapResultRow":
        r = cls()
        r.data["QryContigID"]  = int(qid)
        r.data["RefcontigID1"] = cls.rid(A)
        r.data["RefcontigID2"] = cls.rid(B)

        q_start, q_end, r_start, r_end = cls.span_both_axes(A, B)
        r.data["QryStartPos"] = q_start
        r.data["QryEndPos"]   = q_end
        r.data["RefStartPos"] = r_start
        r.data["RefEndPos"]   = r_end

        if str(sv_type).startswith("translocation"):
            r.data["SVsize"] = -1.0
        else:
            r.data["SVsize"] = float(q_end - q_start)

        r.data["XmapID1"] = cls.get_xmap_id(A)
        r.data["XmapID2"] = cls.get_xmap_id(B)

        if sv_type in ("duplication", "duplication_split", "duplication_inverted"):
            r.data["Confidence"] = -1.0
        else:
            r.data["Confidence"] = 99.0

        r.data["Type"] = sv_type
        if str(sv_type).startswith("translocation"):
            r.data["Orientation"] = orientation_str
        else:
            r.data["Orientation"] = None

        qsi, qei, rsi, rei = cls.edge_label_indices(A, B)
        r.data["QryStartIdx"] = qsi
        r.data["QryEndIdx"]   = qei
        r.data["RefStartIdx"] = rsi
        r.data["RefEndIdx"]   = rei

        return r



@dataclass
class AlignmentResults:
    referenceFilePath: str
    queryFilePath: str
    rows: List[AlignmentResultRow]

    @staticmethod
    def create(referenceFilePath: str,
               queryFilePath: str,
               rows: List[AlignmentResultRow]):
        return AlignmentResults(
            referenceFilePath,
            queryFilePath,
            AlignmentResults.filterOutSubsequentAlignmentsForSingleQuery(rows))

    @staticmethod
    def filterOutSubsequentAlignmentsForSingleQuery(alignmentResultRows):
        rowsSortedByQueryIdThenByConfidence = \
            sorted(sorted(alignmentResultRows, key=lambda r: r.confidence, reverse=True), key=lambda r: r.queryId)
        rowsWithoutSubsequentAlignmentsForSingleQuery = \
            [next(group) for _, group in itertools.groupby(rowsSortedByQueryIdThenByConfidence, lambda r: r.queryId)]
        return rowsWithoutSubsequentAlignmentsForSingleQuery

    @staticmethod
    def resolve(rows: List[AlignmentResultRow],
                maxDifference: int):
        separate = []
        joined = []
        for _, queries in itertools.groupby(sorted(rows, key=lambda r: r.referenceId),
                                            lambda r: r.referenceId):
            for _, group in itertools.groupby(sorted(list(queries), key=lambda r: r.queryId), lambda r: r.queryId):
                group = list(group)
                if len(group) == 1:
                    separate.append(group[0])
                else:
                    if group[0].check_overlap(group[1], maxDifference):
                        resolved = group[0].resolve(group[1])
                        if resolved:
                            joined.append(resolved)
                        else:
                            separate.extend(group)
                    else:
                        separate.extend(group)
        return joined, separate

    def write_indel_file(self, file_name:str="indels.txt"):
        """Function used to write indels files
    
        :param file_name: Name of output file, defaults to "indels.txt"
        :type file_name: str, optional
        """
        lines_sorted = sorted([indel for row in self.rows for indel in row.indelList()],
                        key = lambda indel: (indel[2], indel[4]))
    
        with open(file_name, "w") as f:
            f.write("#Type \t Length \t RefId \t RefStartPos \t RefEndPos \t RefStartIdx \t RefEndIdx \t QueryId \t QueryStartPos \t QueryEndPos \t QueryStartIdx \t QueryEndIdx \n")
            for line in lines_sorted:
                line = [str(i) for i in line]
                line = "\t".join(line)
                f.write(line + "\n")

    def write_rest_file(self, file_name:str="rests.txt"):
        """Function used to write rests files
    
        :param file_name: Name of output file, defaults to "rests.txt"
        :type file_name: str, optional
        """
        items_sorted = sorted([(row.queryId, row.rests) for row in self.rows])
    
        with open(file_name, "w") as f:
            for item in items_sorted:
                rests = item[1]
                rest_maps = 'None' if rests is None else '\t'.join(map(str, rests))
                f.write("{}\t{}\n".format(item[0], rest_maps))







    def debugPrintAlignmentGroup(self, qid: int, group_rows: list) -> None:
        print(f"[SMAP GROUP] qid={qid}  n_alignments={len(group_rows)}")
        for i, r in enumerate(group_rows, 1):
            ori = '-' if getattr(r, 'reverseStrand', False) else '+'
            print(
                f"  #{i}: ref={r.referenceId}  "
                f"qpos=({r.queryStartPosition:.1f},{r.queryEndPosition:.1f}) "
                f"rpos=({r.referenceStartPosition:.1f},{r.referenceEndPosition:.1f}) "
                f"ori={ori}  conf={r.confidence:.3f}"
            )





    def effectiveQStart(self, row) -> float:
        """
        Efektywny 'początek' na query niezależnie od orientacji:
        bierzemy minimum z (queryStartPosition, queryEndPosition).
        Dzięki temu grupy trafiają do SmapResultRow w kolejności rosnącej
        wzdłuż molekuły.
        """
        qs = float(row.queryStartPosition)
        qe = float(row.queryEndPosition)
        return qs if qs <= qe else qe

    def groupAlignmentsByQuery(self) -> Dict[int, List]:
        """
        Zgrupuj wyrównania po molekule (queryId) i posortuj w obrębie grupy
        po faktycznym początku na query (min(qStart, qEnd)), a następnie
        po refId i RefStart (dla stabilności).
        """
        groups: Dict[int, List] = defaultdict(list)
        for r in self.rows:
            groups[r.queryId].append(r)

        for qid, lst in groups.items():
            lst.sort(key=lambda r: (
                self.effectiveQStart(r),
                int(r.referenceId),
                float(r.referenceStartPosition),
            ))
        return groups

    def to_smap_rows(self) -> List["SmapResultRow"]:
        out: List[SmapResultRow] = []
        groups = self.groupAlignmentsByQuery()

        build_group = getattr(SmapResultRow, "from_alignment_group", None)
        has_factory = callable(build_group)

        for qid, group_rows in groups.items():
            self.debugPrintAlignmentGroup(qid, group_rows)

            if has_factory:
                smap_list = build_group(qid, group_rows)

                if smap_list is None:
                    for ar in group_rows:
                        sm = SmapResultRow.from_alignment_row(ar)
                        out.append(sm)
                elif len(smap_list) > 0:
                    out.extend(smap_list)
                else:
                    pass

            else:
                for ar in group_rows:
                    sm = SmapResultRow.from_alignment_row(ar)
                    out.append(sm)

        return out


class AlignmentResultRow(BenchmarkAlignment):
    @staticmethod
    def create(segmentsWithoutConflicts: AlignmentSegmentsWithResolvedConflicts,
               queryId: int,
               referenceId: int,
               queryLength: int,
               referenceLength: int,
               reverseStrand: bool):
        segments = segmentsWithoutConflicts.segments
        alignedPairs = sorted(p for s in segments for p in s.positions if isinstance(p, AlignedPair))
        firstPair = alignedPairs[0] if alignedPairs else AlignedPair.null
        lastPair = alignedPairs[-1] if alignedPairs else AlignedPair.null
        queryStartPosition = (firstPair if not reverseStrand else lastPair).query.position
        queryEndPosition = (lastPair if not reverseStrand else firstPair).query.position
        referenceStartPosition = firstPair.reference.position
        referenceEndPosition = lastPair.reference.position
        confidence = sum(s.segmentScore for s in segments)
        return AlignmentResultRow(segments, queryId, referenceId, queryLength, referenceLength, queryStartPosition,
                                  queryEndPosition, referenceStartPosition, referenceEndPosition, reverseStrand,
                                  confidence)

    def __init__(self,
                 segments: List[AlignmentSegment],
                 queryId: int = 1,
                 referenceId: int = 1,
                 queryLength: int = 1,
                 referenceLength: int = 1,
                 queryStartPosition: int = 0,
                 queryEndPosition: int = 0,
                 referenceStartPosition: int = 0,
                 referenceEndPosition: int = 0,
                 reverseStrand: bool = False,
                 confidence: float = 0.,
                 alignedRest: bool = False,
                 indels: List = [],
                 rests: List | None = None):

        self.queryId = queryId
        self.referenceId = referenceId
        self.queryStartPosition = queryStartPosition
        self.queryEndPosition = queryEndPosition
        self.referenceStartPosition = referenceStartPosition
        self.referenceEndPosition = referenceEndPosition
        self.reverseStrand = reverseStrand
        self.confidence = confidence
        self.queryLength = queryLength
        self.referenceLength = referenceLength
        self.segments = segments
        self.alignedRest = alignedRest
        self.indels = indels
        self.rests = rests

    @property
    def positions(self):
        return [position for segment in self.segments for position in segment.positions]

    @property
    def alignedPairs(self) -> List[AlignedPair]:
        return [p for p in self.positions if isinstance(p, AlignedPair)]

    @property
    def notAlignedPositions(self) -> List[NotAlignedPosition]:
        return [p for p in self.positions if isinstance(p, NotAlignedPosition)]

    @property
    def cigarString(self):
        if not self.alignedPairs:
            return ""
        hitEnums = list(self.__getHitEnums())
        return "".join(self.__aggregateHitEnums(hitEnums))

    def __getHitEnums(self):
        pairs = list(self.__removeDuplicateQueryPositionsPreservingLastOne(self.alignedPairs))
        pairsIterator = iter(pairs)
        currentPair: AlignedPair = next(pairsIterator)
        previousQuery = currentPair.query.siteId
        for referenceIndex in range(pairs[0].reference.siteId,
                                    pairs[-1].reference.siteId + 1):
            queryIncrement = abs(currentPair.query.siteId - previousQuery)
            if queryIncrement > 1:
                for _ in range(1, queryIncrement):
                    yield HitEnum.INSERTION
                previousQuery = currentPair.query.siteId
            if currentPair.reference.siteId == referenceIndex:
                previousQuery = currentPair.query.siteId
                currentPair = next(pairsIterator, None)
                yield HitEnum.MATCH
            elif currentPair.reference.siteId > referenceIndex:
                yield HitEnum.DELETION

    @staticmethod
    def __removeDuplicateQueryPositionsPreservingLastOne(pairs: List[AlignedPair]):
        for _, ambiguousPairs in itertools.groupby(pairs, lambda pair: pair.query.siteId):
            *_, lastPair = ambiguousPairs
            yield lastPair

    @staticmethod
    def __aggregateHitEnums(hits: List[HitEnum]):
        hit = None
        count = 1
        previousHit: HitEnum = hits[0]
        for hit in hits[1:]:
            if hit == previousHit:
                count += 1
            else:
                yield AlignmentResultRow.__hitToString(count, previousHit)
                previousHit = hit
                count = 1
        if hit:
            yield AlignmentResultRow.__hitToString(count, hit)

    @staticmethod
    def __hitToString(count, hit):
        x = f"{count}{hit.value}"
        return x

    def getUnalignedFragments(self, queries: List[OpticalMap]) -> List[OpticalMap]:
        query = next((opticMap for opticMap in queries if opticMap.moleculeId == self.queryId), None)
        if query is None or not self.alignedPairs:
            return []

        start_pos = min(self.queryStartPosition, self.queryEndPosition)
        end_pos = max(self.queryStartPosition, self.queryEndPosition)

        try:
            start_idx = query.positions.index(start_pos)
            end_idx = query.positions.index(end_pos)
        except ValueError:
            return []

        fragments: List[OpticalMap] = []

        if start_idx > 0:
            leading_positions = query.positions[:start_idx]
            fragments.append(OpticalMap(self.queryId, self.queryLength, leading_positions, shift=0))

        if end_idx < len(query.positions) - 1:
            trailing_positions = query.positions[end_idx + 1:]
            fragments.append(
                OpticalMap(
                    self.queryId,
                    self.queryLength,
                    trailing_positions,
                    shift=len(query.positions) - len(trailing_positions)
                )
            )

        return fragments

    def setAlignedRest(self, alignedRest: bool):
        self.alignedRest = alignedRest
        return self

    def check_overlap(self, alignedRest: AlignmentResultRow, maxDifference: int) -> bool:
        """Function used to identify overlapping alignments of the same query

        :param alignedRest: Other alignment of the same query
        :type alignedRest: AlignmentResultRow
        :param maxDifference: Maximum difference between reference positions of the
        alignments if they are to be joint
        :type maxDifference: int
        :return: Whether those two alignments should be joined
        :rtype: bool
        """
        if self.orientation == alignedRest.orientation and self.referenceId == alignedRest.referenceId:
            diff = abs(max(self.referenceStartPosition, alignedRest.referenceStartPosition) - \
                       min(self.referenceEndPosition, alignedRest.referenceEndPosition))
            if diff <= maxDifference:
                return True
        return False

    def resolve(self, alignedRest: AlignmentResultRow) -> AlignmentResultRow | None:
        """Function used to resolve conflicts between two overlapping alignments of the same query

        :param alignedRest: Other alignment of the same query
        :type alignedRest: AlignmentResultRow
        :return: Joint alignment with resolved conflicts
        :rtype: AlignmentResultRow
        """
        if self.alignedPairs[0].reference.position < alignedRest.alignedPairs[0].reference.position:
            pair = self.segments[0].checkForConflicts(alignedRest.segments[0])
        else:
            pair = alignedRest.segments[0].checkForConflicts(self.segments[0])

        if resolution := pair.resolveConflict():
            seg1, seg2 = resolution
            notEmptySegments = [s for s in [seg1, seg2] if s != AlignmentSegment.empty]
            return AlignmentResultRow.create(AlignmentSegmentsWithResolvedConflicts(notEmptySegments),
                                             self.queryId, self.referenceId, self.queryLength, self.referenceLength,
                                             self.reverseStrand)
        return

    def indelList(self):
        indels = []
        for leftAP, rightAP, diff in self.indels:
            refStartPos, queryStartPos, refStartIdx, queryStartIdx = leftAP
            refEndPos, queryEndPos, refEndIdx, queryEndIdx = rightAP
            # diff = abs(refEndPos-refStartPos) - abs(queryEndPos-queryStartPos)
            indelType = 'insertion' if diff<0 else 'deletion'
            indels.append([indelType, diff, 
                        self.referenceId, refStartPos, refEndPos, refStartIdx, refEndIdx, 
                        self.queryId,  queryStartPos, queryEndPos,  queryStartIdx, queryEndIdx])
        return indels


    # def indelList(self):
    #     if len(self.segments) <= 1:
    #         return []
    #     indels = []
    #     currAPs = [ap for ap in self.segments[0].positions if isinstance(ap, AlignedPair)]
    #     leftAP = currAPs[-1]
    #     for currSegment in self.segments[1:]:
    #         currAPs = [ap for ap in currSegment.positions if isinstance(ap, AlignedPair)]
    #         rightAP = currAPs[0]
    #         refStartPos = leftAP.reference.position
    #         refEndPos = rightAP.reference.position
    #         queryStartPos = leftAP.query.position
    #         queryEndPos = rightAP.query.position
    #         refStartIdx = leftAP.reference.siteId
    #         refEndIdx = rightAP.reference.siteId
    #         queryStartIdx = leftAP.query.siteId
    #         queryEndIdx = rightAP.query.siteId
    #         diff = abs(refEndPos-refStartPos) - abs(queryEndPos-queryStartPos)
    #         indelType = 'insertion' if diff<0 else 'deletion'
    #         indels.append([indelType, diff, 
    #                    self.referenceId, refStartPos, refEndPos, refStartIdx, refEndIdx, 
    #                    self.queryId,  queryStartPos, queryEndPos,  queryStartIdx, queryEndIdx])
    #         leftAP = currAPs[-1]
    #     return indels



