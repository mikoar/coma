import itertools
import pandas
import re
from maps.reference_optical_map import ReferenceOpticalMap
from maps.sequence_generator import SequenceGenerator
from maps.optical_map import OpticalMap, PositionRange


def readOpticalMap(sdataFilePath: str, resolution: int, blurRadius: int, moleculeId: int):
    return __readSdata(sdataFilePath, SequenceGenerator(resolution, blurRadius), moleculeId)


def readReference(cmapFilePath: str, resolution: int, blurRadius: int):
    return __readReferenceCmap(cmapFilePath, SequenceGenerator(resolution, blurRadius))


# def __readOpticalMap(filePath, builder):
#     extension = os.path.splitext(filePath)[1][1:]

#     if extension == "cmap":
#         return __readCmap(filePath, moleculeId),
#     elif extension == "sdata":
#         return __readSdata(filePath, moleculeId)
#     else:
#         raise ValueError(extension)


def __readReferenceCmap(filePath: str, sequenceGenerator: SequenceGenerator):
    maps = pandas.read_csv(
        filePath, comment="#", delimiter="\t", names=__getCmapColumnNames(filePath))

    # if moleculeId:
    #     maps = maps[maps["CMapId"] == moleculeId]

    positions = [int(position) for _, position in maps["Position"].sort_values().iteritems()]
    sequence = sequenceGenerator.positionsToSequence(positions)
    return ReferenceOpticalMap(sequence)


def __getCmapColumnNames(filePath):
    with open(filePath) as file:
        gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
        header_line = list(itertools.islice(gen, 1))[0].strip()
        names = re.split('\s+', header_line)[1:]
    return names


def __readSdata(filePath, sequenceGenerator: SequenceGenerator,  moleculeId):
    if not moleculeId:
        raise ValueError(moleculeId)

    names = ["Fragment ID",
             "Reference",
             "Strand",
             "Start",
             "Stop",
             "SimuInfoDetail",
             "Size",
             "TotalSegments",
             "SegmentDetail"]

    map = pandas.read_csv(filePath, comment="#", delimiter="\t", names=names, index_col="Fragment ID").loc[moleculeId]
    referenceCoordinates = PositionRange(int(map["Start"]), int(map["Stop"]))
    segments = [int(segment) for segment in str(map["SegmentDetail"]).split(";")]
    sequence = sequenceGenerator.segmentsToSequence(segments)

    return OpticalMap(moleculeId, map["Strand"], sequence, referenceCoordinates)
