import sys

sys.path.insert(0, "../")
from collections import defaultdict
from src.parsers.cmap_reader import CmapReader
from src.parsers.xmap_reader import XmapReader
import pandas as pd


def read_all_files(reference_file:str, alignment_file:str, query_file:str) -> tuple[dict, dict, dict]:
    """Function used to read all files used for alignment

    :param reference_file: name of reference file
    :type reference_file: str
    :param alignment_file: name of alignment file
    :type alignment_file: str
    :param query_file: name of query file
    :type query_file: str
    :return: dictionaries with parsed molecules and alignments
    :rtype: tuple[dict, dict, dict]
    """
    alignments = XmapReader().readAlignments(open(alignment_file))
    reference = CmapReader().readReferences(open(reference_file))
    query = CmapReader().readReferences(open(query_file))

    alignment_queryId = {}
    alignment_referenceId = {}
    alignment_queryId = defaultdict(lambda: [], alignment_queryId)
    alignment_referenceId = defaultdict(lambda: [], alignment_referenceId)
    reference_dict = {}
    query_dict = {}

    for alignment in alignments:
        alignment_referenceId[alignment.referenceId].append(alignment)
        alignment_queryId[alignment.queryId].append(alignment)
    for ref in reference:
        if ref:
            reference_dict[ref.moleculeId] = ref
    for quer in query:
        if quer:
            query_dict[quer.moleculeId] = quer
    return reference_dict, query_dict, alignment_referenceId


def read_alignments_file(alignment_file:str) -> dict:
    """Function used to read alignment file and returned parsed dictionary

    :param alignment_file: Name of aligned XMAP file
    :type alignment_file: str
    :return: Dictionary with reference ids as keys and alignments as values
    :rtype: dict
    """
    alignment_referenceId = {}
    alignment_referenceId = defaultdict(lambda: [], alignment_referenceId)
    alignments = XmapReader().readAlignments(open(alignment_file))
    for alignment in alignments:
        alignment_referenceId[alignment.referenceId].append(alignment)
    return  alignment_referenceId


def read_segments_file(seg_catch_file: str) -> pd.DataFrame:
    """Function used to read csv file with molecule segments

    :param seg_catch_file: Name of the COMA file with created segments
    :type seg_catch_file: str
    :return: DataFrame with segments created during COMA workflow
    :rtype: pd.DataFrame
    """
    df = pd.read_csv(seg_catch_file, sep=";")
    df = df.sort_values(by="queryId")
    small_df = df[~df["segmentNb"].isin([0,1])]
    small_df = small_df[ small_df.groupby(['queryId'])['alignmentConfidence'].transform(max) == small_df['alignmentConfidence']]
    return small_df
        