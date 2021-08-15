# %%
import pandas as pd
from matplotlib import pyplot as plt
from pandas.core.groupby.groupby import GroupBy

from src.cmap_reader import AlignmentReader, BionanoFileReader, CmapReader
from src.plot import plotCorrelation
from src.sequence_generator import SequenceGenerator

alignmentsFile = "data/NA12878_BSPQI/EXP_REFINEFINAL1.xmap"
referenceFile = "data/NA12878_BSPQI/hg19_NT.BSPQI_0kb_0labels.cmap"
queryFile = "data/NA12878_BSPQI/EXP_REFINEFINAL1.cmap"


def plot(alignmentIds, resolution=128, blur=4):
    alignmentReader = AlignmentReader()
    alignments = alignmentReader.readAlignments(alignmentsFile, alignmentIds)
    queryIds = [a.queryId for a in alignments]
    referenceIds = [a.referenceId for a in alignments]
    sequenceGenerator = SequenceGenerator(resolution, blur)
    sequenceReader = CmapReader(sequenceGenerator)
    queries = sequenceReader.readQueries(queryFile, queryIds)
    references = sequenceReader.readReferences(referenceFile, referenceIds)
    for alignment in alignments:
        query = next(q for q in queries if q.moleculeId == alignment.queryId)
        reference = next(r for r in references if r.moleculeId == alignment.referenceId)
        result = query.correlate(reference.sequence, reverseStrand=alignment.reverseStrand)
        fig = plotCorrelation(result, resolution, (alignment.expectedQueryMoleculeStart, alignment.expectedQueryMoleculeEnd))
        fig.suptitle(f'Alignment {alignment.id}')
        fig.savefig(f"output_heatmap/not_mapped_molecules/alignmentId_{alignment.id}.svg", bbox_inches='tight', pad_inches=0)


def validityRatio(isValidColumn: pd.Series):
    return isValidColumn.sum() / isValidColumn.size


def getMappedRatio(group: GroupBy):
    return pd.Series({'mappedRatio': validityRatio(group['isValid'])})


def getAlignmentLengthToQueryLength(group: GroupBy):
    return pd.Series({'alignmentLengthToQueryLength':
                      abs(group['QryEndPos'].iloc[0] - group['QryStartPos'].iloc[0]) / group['QryLen'].iloc[0],
                      'length':  group['QryLen'].iloc[0]})


def getAlignmentLength(group: GroupBy):
    return pd.Series({'alignmentLength': min(
                      abs(group['QryEndPos'].iloc[0] - group['QryStartPos'].iloc[0]),
                      abs(group['RefEndPos'].iloc[0] - group['RefStartPos'].iloc[0]))})


def getConfidence(group: GroupBy):
    return pd.Series({'confidence': group['confidence'].iloc[0]})


results = pd.read_csv(
    "data/result_peakWithinAlignmentSizeFromCenter_count_1703_res_64,128,256,512,1024_blur_0,2,4,8,16.csv").set_index(['resolution', 'blur'])
bestResBlurPairs = [(64, 4), (64, 8), (128, 2), (128, 4), (256, 2)]
results = results[results.index.isin(bestResBlurPairs)]

# %%
groupedByAlignment = results.groupby('alignmentId')
notAlignedAnywhere = groupedByAlignment.filter(lambda group: validityRatio(group['isValid']) <= 0.).reset_index()
notAlignedAlignmentIds = notAlignedAnywhere['alignmentId'].unique().tolist()
plot(notAlignedAlignmentIds)

# %%


groupedByAlignment = results.groupby('alignmentId')
mappedRatio = groupedByAlignment.apply(getMappedRatio)
mappedRatio.to_csv('output_heatmap/not_mapped_molecules/alignmentMappedRatio.csv')

# %%


alignments = BionanoFileReader().readFile(alignmentsFile, ["XmapEntryID", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos", "QryLen"])
resultsWithLengths = results.reset_index().join(alignments.set_index('XmapEntryID'), on='alignmentId', how='left')
groupedByAlignment = resultsWithLengths.groupby('alignmentId')

mappedRatio: pd.DataFrame = groupedByAlignment.apply(getMappedRatio)
alignmentLengthToQueryLength: pd.DataFrame = groupedByAlignment.apply(getAlignmentLengthToQueryLength)

queryAlignmentSpanVsMappingRatio = pd.concat([mappedRatio, alignmentLengthToQueryLength], axis=1)
queryAlignmentSpanVsMappingRatio[['mappedRatio', 'alignmentLengthToQueryLength']].plot(
    x='mappedRatio', y='alignmentLengthToQueryLength', kind='scatter')
plt.savefig('output_heatmap/not_mapped_molecules/mapped_ratio_vs_alignment_length_to_query_length.svg')


# %%
alignmentLength: pd.DataFrame = groupedByAlignment.apply(getAlignmentLength)
alignmentLengthVsMappingRatio = pd.concat([mappedRatio, alignmentLength], axis=1)
alignmentLengthVsMappingRatio[['mappedRatio', 'alignmentLength']].plot(
    x='mappedRatio', y='alignmentLength', kind='scatter')
# %%
confidence: pd.DataFrame = groupedByAlignment.apply(getConfidence)
confidenceVsMappingRatio = pd.concat([mappedRatio, confidence], axis=1)
confidenceVsMappingRatio[['mappedRatio', 'confidence']].plot(
    x='mappedRatio', y='confidence', kind='scatter')
# %%
qualityHist = plt.hist(results[results.index == (128, 4)].score, 100)

# przeanalizować niezmapowane contigi, + porównnać długość alignmentu/długość query, scharakteryzować dlaczego się nie mapują
# znaleźć contigi, które nigdzie się nie mapują (przy żadnych parametrach)
# wziąć 8 kafelków >98%,
# podzielić contigi na grupy - te co się lepiej/gorzej mapują, jaka część contigu jest zawarta w alignmencie
# odrzucić w umotywowany sposób słabe alignmenty, zrobić heatmapę bez nich
# wykres rozkład stosunku długości alignmentu do długości zapytania i jak to się ma do wykrywalności
# lista alignmentId -> % zmapowane
# %%

# pełne heatmapy dla 2 nwych funkcji waidacji
# wybrać najlepsze kafelki, tylko je do scatter plotów, y - skala log
# może plot 3d dla lepszego  odseparowania niealignujących się cząsteczek
# histogram tylko dla wybranych kafelków
# powtózenie dla molekuł, nie contigów (NA12878_BSPQI_pipeline_results/output/contigs/alignmolvref/merge/)
