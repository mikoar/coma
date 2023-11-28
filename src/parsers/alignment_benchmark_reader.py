import itertools
from typing import TextIO, Iterable, List

from src.diagnostic.benchmark_alignment import BenchmarkAlignment
from src.parsers.simulation_data_as_xmap_reader import SimulationDataAsXmapReader
from src.parsers.xmap_reader import XmapReader


class AlignmentBenchmarkReader:
    def __init__(self, xmapReader: XmapReader, simulationReader: SimulationDataAsXmapReader):
        self.xmapReader = xmapReader
        self.simulationReader = simulationReader

    def read(self, file: TextIO, queryIds: Iterable[int] = None) -> List[BenchmarkAlignment]:
        headers = "\n".join(itertools.islice(itertools.takewhile(lambda line: line.startswith("#"), file), 10))
        file.seek(0, 0)
        if "XMAP" in headers:
            return self.xmapReader.readAlignments(file, queryIds=queryIds)
        if "SimuInfoDetail" in headers:
            return self.simulationReader.readAlignments(file, queryIds=queryIds)
        else:
            raise Exception(f"File {file.name} is in unknown format. Either XMAP or SDATA formats are supported.")
