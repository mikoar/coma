# %%
import pandas as pd

from cmap_reader import AlignmentReader, CmapReader
from plot import plotCorrelation
from sequence_generator import SequenceGenerator

alignmentsFile = "../data/EXP_REFINEFINAL1.xmap"
referenceFile = "../data/hg19_NT.BSPQI_0kb_0labels.cmap"
queryFile = "../data/EXP_REFINEFINAL1.cmap"


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
        fig = plotCorrelation(result, resolution, (alignment.expectedQueryStart, alignment.expectedQueryEnd))
        fig.savefig(f"../output_stats/not_mapped_molecules/alignmentId_{alignment.id}.svg", bbox_inches='tight', pad_inches=0)


def maxValidityRatio(isValidColumn: pd.Series, maxValidityRatio: float):
    return isValidColumn.sum() / isValidColumn.size <= maxValidityRatio


results = pd.read_csv("../output_stats/result_count_1000_res_32,48,64,128,256,512_blur_0,2,4,8,16.csv").set_index(['resolution', 'blur', 'alignmentId'])

# %%
groupedByAlignment = results.groupby('alignmentId')
notAlignedAnywhere = groupedByAlignment.filter(lambda group: maxValidityRatio(group['isValid'], 0.)).reset_index()
notAlignedAlignmentIds = notAlignedAnywhere['alignmentId'].unique().tolist()
plot(notAlignedAlignmentIds)

# przeanalizować niezmapowane contigi, + porównnać długość alignmentu/długość query, scharakteryzować dlaczego się nie mapują
# znaleźć contigi, które nigdzie się nie mapują (przy żadnych parametrach)
# wziąć 8 kafelków >98%,
# podzielić contigi na grupy - te co się lepiej/gorzej mapują, jaka część contigu jest zawarta w alignmencie
# odrzucić w umotywowany sposób słabe alignmenty, zrobić heatmapę bez nich
# wykres rozkład stosunku długości alignmentu do długości zapytania i jak to się ma do wykrywalności
# lista alignmentId -> % zmapowane
# %%
