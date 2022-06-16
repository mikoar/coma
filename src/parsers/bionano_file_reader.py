from __future__ import annotations

import itertools
import re
from typing import List, TextIO

import pandas
from pandas import DataFrame


class BionanoFileReader:
    def readFile(self, file: TextIO, columns: List[str]) -> DataFrame:
        return pandas.read_csv(
            file,
            comment="#",
            delimiter="\t",
            names=self.__getColumnNames(file),
            usecols=columns)

    @staticmethod
    def __getColumnNames(file: TextIO):
        commentLines = itertools.dropwhile(lambda line: not line.startswith('#h'), file)
        header_line = list(itertools.islice(commentLines, 1))[0].strip()
        names = re.split(r'\s+', header_line)[1:]
        return names
