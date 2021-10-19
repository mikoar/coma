import itertools
import re
from typing import List

import pandas
from pandas import DataFrame

from src.alignment.alignment_results import AlignmentResults


class BionanoFileReader:
    def readFile(self, filePath, columns: List[str]) -> DataFrame:
        return pandas.read_csv(
            filePath,
            comment="#",
            delimiter="\t",
            names=self.__getColumnNames(filePath),
            usecols=columns)

    # def writeFile(self, filePath, alignmentResults: AlignmentResults):
    #     dataFrame = DataFrame.
    #     pandas.to_csv(
    #         filePath,
    #         comment="#",
    #         delimiter="\t",
    #         names=self.__getColumnNames(filePath),
    #         usecols=columns)

    def __getColumnNames(self, filePath):
        with open(filePath) as file:
            gen = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
            header_line = list(itertools.islice(gen, 1))[0].strip()
            names = re.split(r'\s+', header_line)[1:]
        return names
