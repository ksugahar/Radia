"""Pipe framing, native-import transport lifecycle, and fail-loud protection."""

import ctypes
from ctypes import wintypes
import io
import sys
import threading
from types import SimpleNamespace

import pytest

from radia_mcp.common import windows_stdio as transport


def _set(pointer, value):
    ctypes.cast(pointer, ctypes.POINTER(wintypes.DWORD)).contents.value = value


class Pipe:
    def __init__(self, data=b"", *, chunk=3, idle=False):
        self.data = bytearray(data)
        self.chunk = chunk
        self.idle = idle
        self.waiting = threading.Event()
        self.requests = []

    def PeekNamedPipe(self, handle, unused, size, read, available, left):
        if not self.data and not self.idle:
            return False
        _set(available, min(len(self.data), self.chunk))
        if not self.data:
            self.waiting.set()
        return True

    def ReadFile(self, handle, target, count, received, overlapped):
        assert 0 < count <= min(len(self.data), self.chunk)
        self.requests.append(count)
        ctypes.memmove(target, bytes(self.data[:count]), count)
        del self.data[:count]
        _set(received, count)
        return True


@pytest.mark.parametrize("chunk", [3, 65536])
@pytest.mark.parametrize("eof_error", [109, 232, 233])
def test_fragmented_utf8_lines_and_eof_preserve_bytes(monkeypatch, chunk, eof_error):
    monkeypatch.setattr(ctypes, "get_last_error", lambda: eof_error, raising=False)
    payload = ('{"text":"場の確認"}\n' + '{"data":"' + "x" * 20000 + '"}\n').encode()
    pipe = Pipe(payload, chunk=chunk)
    raw = transport._PipeReader(123, pipe, 0)
    with io.TextIOWrapper(io.BufferedReader(raw), encoding="utf-8") as stream:
        assert "".join(stream).encode() == payload
    assert max(pipe.requests) <= min(pipe.chunk, io.DEFAULT_BUFFER_SIZE)
    assert raw.closed


@pytest.mark.parametrize("error", [5, 6])
def test_unexpected_pipe_errors_are_not_eof(monkeypatch, error):
    monkeypatch.setattr(ctypes, "get_last_error", lambda: error, raising=False)
    raw = transport._PipeReader(123, Pipe(), 0)
    with pytest.raises(OSError) as failure:
        raw.readinto(bytearray(4))
    assert failure.value.errno == error
    raw.close()


def test_scoped_failure_restores_stdin_and_releases_pending_reader(monkeypatch):
    previous = sys.stdin
    pipe = Pipe(idle=True)
    raw = transport._PipeReader(123, pipe, 0)
    monkeypatch.setattr(transport, "_reader_for_stdin", lambda stream: raw)
    result = []
    worker = None
    with pytest.raises(RuntimeError, match="tool failure"):
        with transport._pipe_stdin():
            stream = sys.stdin
            worker = threading.Thread(target=lambda: result.append(stream.readline()))
            worker.start()
            assert pipe.waiting.wait(1)
            raise RuntimeError("tool failure")
    worker.join(1)
    assert not worker.is_alive()
    assert result == [""]
    assert sys.stdin is previous and not previous.closed
    raw.close()  # Idempotent; never closes the borrowed OS handle/descriptor.


def test_non_pipe_uses_original_sdk_transport(monkeypatch):
    previous = sys.stdin
    monkeypatch.setattr(transport, "_reader_for_stdin", lambda stream: None)
    calls = []
    server = SimpleNamespace(run=lambda **options: calls.append((options, sys.stdin)))
    transport.run_stdio(server)
    assert calls == [({"transport": "stdio"}, previous)]


def test_closed_reader_rejects_reads():
    raw = transport._PipeReader(123, Pipe(), 0)
    raw.close()
    with pytest.raises(ValueError, match="closed"):
        raw.readinto(bytearray(1))
