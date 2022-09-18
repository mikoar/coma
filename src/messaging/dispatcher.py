from typing import List

from src.messaging.message_handler import MessageHandler
from src.messaging.messages import Message


class Dispatcher:
    def __init__(self, handlers: List[MessageHandler] = None):
        self.__handlers = handlers or []

    def addHandler(self, handler: MessageHandler):
        self.__handlers.append(handler)

    def dispatch(self, message: Message):
        for handler in self.__handlers:
            if handler.canHandle(message):
                handler.handle(message)
