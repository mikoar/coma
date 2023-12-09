from typing import List

from src.extensions.extension import Extension
from src.extensions.messages import Message


class Dispatcher:
    def __init__(self, extensions: List[Extension] = None):
        self.__extensions = extensions or []

    def addExtension(self, extension: Extension):
        self.__extensions.append(extension)

    def dispatch(self, message: Message):
        for extension in self.__extensions:
            if extension.canHandle(message):
                extension.handle(message)
