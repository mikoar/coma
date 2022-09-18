from abc import ABC

from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.optical_map import InitialAlignment, CorrelationResult


class Message(ABC):
    def type(self):
        return type(self)


class InitialAlignmentMessage(Message):
    def __init__(self, data: InitialAlignment):
        self.data = data


class CorrelationResultMessage(Message):
    def __init__(self, initialAlignment: InitialAlignment, refinedAlignment: CorrelationResult):
        self.initialAlignment = initialAlignment
        self.refinedAlignment = refinedAlignment


class AlignmentResultRowMessage(Message):
    def __init__(self, data: AlignmentResultRow):
        self.data = data
