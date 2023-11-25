from __future__ import annotations

import itertools
import re
from typing import List, TextIO

import pandas
from pandas import DataFrame


class BionanoFileReader:
    def __init__(self, headersLinePrefix: str = '#h'):
        self.headersLinePrefix = headersLinePrefix

    def readFile(self, file: TextIO, columns: List[str]) -> DataFrame:
        return pandas.read_csv(
            file,
            comment="#",
            delimiter="\t",
            names=self.__getColumnNames(file),
            usecols=columns)

    def __getColumnNames(self, file: TextIO):
        commentLines = itertools.dropwhile(lambda line: not line.startswith(self.headersLinePrefix), file)
        header_line = list(itertools.islice(commentLines, 1))[0].strip()
        names = re.split(r'\s+', header_line)[1:]
        return names
