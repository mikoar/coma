from abc import abstractmethod
from typing import Type

from src.messaging.messages import Message


class MessageHandler:
    @property
    @abstractmethod
    def messageType(self) -> Type[Message]:
        pass

    def canHandle(self, message: Message):
        return message.type() == self.messageType

    @abstractmethod
    def handle(self, message: Message):
        pass
