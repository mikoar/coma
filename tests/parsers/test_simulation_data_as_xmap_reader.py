from typing import Dict, List
from unittest.mock import Mock

import pytest
from pandas import DataFrame

from src.diagnostic.benchmark_alignment import BenchmarkAlignedPair
from src.parsers.bionano_file_reader import BionanoFileReader
from src.parsers.simulation_data_as_xmap_reader import SimulationDataAsXmapReader


def test_parsesMap():
    sdataReader = __getSut({
        "ID": ["51", "10"],
        "Reference": ["1", "2"],
        "Strand": ["+", "-"],
        "Start": ["4147184", "105011676"],
        "Stop": ["4347183", "105211675"],
        "SimuInfoDetail": ["1:657;1:658;1:659;1:660;FP;1:661;FP;1:662;1:663;1:664;1:665;FP;1:666;1:667;1:669;1:670;1:671;1:672;1:673",
                           "2:14649,2:14650;FP;2:14648;2:14647;2:14646;2:14645;2:14644;FP;2:14643;2:14642;2:14641"],
        "Size": ["199860", "200048"],
        "TotalSegments": ["20", "12"],
        "SegmentDetail": ["9117;10816;25859;12554;20475;2595;2022;9595;4160;23252;3610;5944;979;2641;25903;6710;2719;4432;8454;18004",
                          "11664;881;45688;4466;12162;14678;16875;25870;2175;5388;3855;56335"]})

    xmaps = sdataReader.readAlignments(Mock())

    assert len(xmaps) == 2
    xmap1 = xmaps[0]
    assert xmap1.queryId == 51
    assert xmap1.referenceId == 1
    assert xmap1.orientation == "+"
    assert xmap1.queryLength == 199860
    assert xmap1.referenceStartPosition == 4147184
    assert xmap1.referenceEndPosition == 4347183
    assert xmap1.alignedPairs == [BenchmarkAlignedPair.create(str(r), str(q)) for r, q in [
        (657, 1), (658, 2), (659, 3), (660, 4), (661, 6), (662, 8), (663, 9), (664, 10), (665, 11), (666, 13), (667, 14), (669, 15), (670, 16), (671, 17),
        (672, 18), (673, 19)
    ]]

    xmap2 = xmaps[1]
    assert xmap2.queryId == 10
    assert xmap2.referenceId == 2
    assert xmap2.orientation == "-"
    assert xmap2.queryLength == 200048
    assert xmap2.referenceStartPosition == 105011676
    assert xmap2.referenceEndPosition == 105211675
    assert xmap2.alignedPairs == [BenchmarkAlignedPair.create(str(r), str(q)) for r, q in [
        (14641, 11), (14642, 10), (14643, 9), (14644, 7), (14645, 6), (14646, 5), (14647, 4), (14648, 3), (14650, 1), (14649, 1)
    ]]


def __getSut(data: Dict[str, List[str]]):
    fileReaderMock: BionanoFileReader = Mock(spec=BionanoFileReader)
    fileReaderMock.readFile = lambda _, __: DataFrame(data=data)
    return SimulationDataAsXmapReader(pairParser=None, reader=fileReaderMock)


if __name__ == '__main__':
    pytest.main(args=[__file__])
