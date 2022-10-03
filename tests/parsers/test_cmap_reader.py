from typing import Dict, List
from unittest.mock import Mock

import pytest
from pandas import DataFrame

from src.parsers.bionano_file_reader import BionanoFileReader
from src.parsers.cmap_reader import CmapReader


@pytest.mark.parametrize("data", [
    ({"CMapId": [1, 1, 1],
      "Position": [100, 200, 300],
      "LabelChannel": [1, 1, 0]}),
    pytest.param(
        {"CMapId": [1, 1, 1, 2],
         "Position": [100, 200, 300, 123],
         "LabelChannel": [1, 1, 0, 0]},
        id="with empty map that should be filtered out")
])
def test_parsesMap(data):
    cmapReader = __getSut(data)

    opticalMaps = cmapReader.readQueries(Mock())

    assert len(opticalMaps) == 1
    opticalMap = opticalMaps[0]
    assert opticalMap.moleculeId == 1
    assert opticalMap.length == 300
    assert opticalMap.positions == [100, 200]


def test_parsesMap_emptyMap():
    data = {
        "CMapId": [1],
        "Position": [100],
        "LabelChannel": [0]}
    cmapReader = __getSut(data)
    opticalMaps = cmapReader.readQueries(Mock())

    assert len(opticalMaps) == 0


def __getSut(data: Dict[str, List[int]]):
    fileReaderMock: BionanoFileReader = Mock(spec=BionanoFileReader)
    fileReaderMock.readFile = lambda _, __: DataFrame(data=data)
    return CmapReader(fileReaderMock)


if __name__ == '__main__':
    pytest.main(args=[__file__])
