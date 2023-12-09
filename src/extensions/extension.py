from abc import abstractmethod
from typing import Type

from src.extensions.messages import Message


class Extension:
    @property
    @abstractmethod
    def messageType(self) -> Type[Message]:
        pass

    def canHandle(self, message: Message):
        return message.type() == self.messageType

    @abstractmethod
    def handle(self, message: Message):
        pass
