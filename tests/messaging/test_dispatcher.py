import pytest

from src.extensions.dispatcher import Dispatcher
from src.extensions.extension import Extension
from src.extensions.messages import Message


class _TestMessage1(Message):
    pass


class _TestMessage2(Message):
    pass


class _TestHandler1(Extension):
    messageType = _TestMessage1
    handledMessages = []

    def handle(self, message: Message):
        self.handledMessages.append(message)


class _TestHandler2(Extension):
    messageType = _TestMessage2
    handledMessages = []

    def handle(self, message: Message):
        self.handledMessages.append(message)


def test_dispatches_messages_to_correct_handlers():
    handler1 = _TestHandler1()
    handler2 = _TestHandler2()
    dispatcher = Dispatcher([handler1, handler2])
    message1 = _TestMessage1()
    message2 = _TestMessage2()

    dispatcher.dispatch(message1)
    dispatcher.dispatch(message2)

    assert len(handler1.handledMessages) == 1
    assert len(handler2.handledMessages) == 1
    assert handler1.handledMessages[0] == message1
    assert handler2.handledMessages[0] == message2


if __name__ == '__main__':
    pytest.main(args=[__file__])
