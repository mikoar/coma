# %%
import pandas as pd
from matplotlib import pyplot as plt, ticker
from pandas.core.groupby.groupby import GroupBy

from src.correlation.cmap_reader import AlignmentReader, CmapReader, BionanoFileReader
from src.correlation.plot import plotCorrelation
from src.correlation.sequence_generator import SequenceGenerator

baseDir = '.local_data/NA12878_BSPQI_pipeline_results/output/contigs/alignmolvref/merge/alignmolvref_contig21'
alignmentsFile = f"{baseDir}.xmap"
referenceFile = f"{baseDir}_r.cmap"
queryFile = f"{baseDir}_q.cmap"


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
        fig = plotCorrelation(result, resolution,
                              (alignment.expectedQueryMoleculeStart, alignment.expectedQueryMoleculeEnd))
        fig.suptitle(f'Alignment {alignment.alignmentId}')
        fig.savefig(f"output_heatmap/not_mapped_molecules/alignmentId_{alignment.alignmentId}.svg", bbox_inches='tight',
                    pad_inches=0)


def validityRatio(isValidColumn: pd.Series):
    return isValidColumn.sum() / isValidColumn.size


def getMappedRatio(group: GroupBy):
    return pd.Series({'mappedRatio': validityRatio(group['isValid'])})


def getAlignmentLengthToQueryLength(group: GroupBy):
    return pd.Series({'alignmentLengthToQueryLength':
                      abs(group['QryEndPos'].iloc[0] - group['QryStartPos'].iloc[0]) / group['QryLen'].iloc[0],
                      'length': group['QryLen'].iloc[0]})


def getAlignmentLength(group: GroupBy):
    return pd.Series({'alignmentLength': min(
        abs(group['QryEndPos'].iloc[0] - group['QryStartPos'].iloc[0]),
        abs(group['RefEndPos'].iloc[0] - group['RefStartPos'].iloc[0]))})


def getConfidence(group: GroupBy):
    return pd.Series({'confidence': group['confidence'].iloc[0]})


title = "contig21_count_9127_res_128,256,512,1024_blur_1,2,3,4"
fileName = f"output_heatmap/result_{title}.csv"
# fileName = "data/result_peakWithinAlignmentSizeFromCenter_count_1703_res_64,128,256,512,1024_blur_0,2,4,8,16.csv"
# fileName = "output_heatmap/result_peakWithinAlignmentSizeUncertainityFromCenterWithFixedMargin_count_1703_res_128,256,512,1024_blur_1,2,3,4.csv"
results = pd.read_csv(fileName).set_index(['resolution', 'blur'])
# bestResBlurPairs = [(64, 4), (64, 8), (128, 2), (128, 4), (256, 2)]
# results = results[results.index.isin(bestResBlurPairs)]

# %%

# %%
# groupedByAlignment = results.groupby('alignmentId')
# notAlignedAnywhere = groupedByAlignment.filter(lambda group: validityRatio(group['isValid']) <= 0.).reset_index()
# notAlignedAlignmentIds = notAlignedAnywhere['alignmentId'].unique().tolist()
# plot(notAlignedAlignmentIds)

# %%


groupedByAlignment = results.groupby('alignmentId')
mappedRatio = groupedByAlignment.apply(getMappedRatio)
# mappedRatio.to_csv('output_heatmap/not_mapped_molecules/alignmentMappedRatio.csv')
#
# # %%
#
#
alignments = BionanoFileReader().readFile(alignmentsFile,
                                          ["XmapEntryID", "QryStartPos", "QryEndPos", "RefStartPos", "RefEndPos",
                                           "QryLen"])
resultsWithLengths = results.reset_index().join(alignments.set_index('XmapEntryID'), on='alignmentId', how='left')
groupedByAlignment = resultsWithLengths.groupby('alignmentId')
#
mappedRatio: pd.DataFrame = groupedByAlignment.apply(getMappedRatio)
alignmentLengthToQueryLength: pd.DataFrame = groupedByAlignment.apply(getAlignmentLengthToQueryLength)
#
queryAlignmentSpanVsMappingRatio = pd.concat([mappedRatio, alignmentLengthToQueryLength], axis=1)
queryAlignmentSpanVsMappingRatio[['mappedRatio', 'alignmentLengthToQueryLength']].plot(
    x='mappedRatio', y='alignmentLengthToQueryLength', kind='scatter')
plt.savefig(f'output_heatmap/plots/queryAlignmentSpanVsMappingRatio_{title}.svg')

# %%
alignmentLength: pd.DataFrame = groupedByAlignment.apply(getAlignmentLength)
alignmentLengthVsMappingRatio = pd.concat([mappedRatio, alignmentLength], axis=1)
ax = alignmentLengthVsMappingRatio[['mappedRatio', 'alignmentLength']].plot(
    x='mappedRatio', y='alignmentLength', kind='scatter', logy=True)
ax.figure.savefig(f'output_heatmap/plots/alignmentLengthVsMappingRatio_{title}.svg')
# %%
confidence: pd.DataFrame = groupedByAlignment.apply(getConfidence)
confidenceVsMappingRatio = pd.concat([mappedRatio, confidence], axis=1)
ax = confidenceVsMappingRatio[['mappedRatio', 'confidence']].plot(
    x='mappedRatio', y='confidence', kind='scatter', logy=True)
ax.figure.savefig(f'output_heatmap/plots/confidenceVsMappingRatio_{title}.svg')
# %%

resolutions, blurs = list(zip(*results.index.drop_duplicates()))
resolutionsCount = len(set(resolutions))
blursCount = len(set(blurs))
fig, axes = plt.subplots(resolutionsCount, blursCount, figsize=(blursCount * 2.5, resolutionsCount * 2.5),
                         constrained_layout=True)
xTicker = ticker.MultipleLocator(0.4)
for resolution, blur, ax in zip(resolutions, blurs, axes.flat):
    ax.hist(results[results.index == (resolution, blur)].score, 100, range=(-0.2, 1))
    ax.xaxis.set_major_locator(xTicker)
    ax.set_title(f"resolution: {resolution}, blur: {blur}")

fig.savefig(f'output_heatmap/plots/hist_{title}.svg')
