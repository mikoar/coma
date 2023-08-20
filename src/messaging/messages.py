from abc import ABC

from src.alignment.alignment_results import AlignmentResultRow
from src.correlation.optical_map import InitialAlignment, CorrelationResult, OpticalMap


class Message(ABC):
    def type(self):
        return type(self)


class InitialAlignmentMessage(Message):
    def __init__(self, data: InitialAlignment):
        self.data = data


class CorrelationResultMessage(Message):
    def __init__(self, initialAlignment: InitialAlignment, refinedAlignment: CorrelationResult, index: int = 0):
        self.initialAlignment = initialAlignment
        self.refinedAlignment = refinedAlignment
        self.index = index


class AlignmentResultRowMessage(Message):
    def __init__(self, reference: OpticalMap, query: OpticalMap, alignment: AlignmentResultRow,
                 correlation: InitialAlignment, index: int = 0):
        self.reference = reference
        self.query = query
        self.alignment = alignment
        self.correlation = correlation
        self.index = index
