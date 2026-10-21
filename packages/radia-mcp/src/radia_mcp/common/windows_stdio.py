"""Keep idle Windows pipe reads from blocking lazy native-module loading.

A pending synchronous stdin read can stall native DLL initialization on Windows.
Read only bytes already available in the pipe. The SDK still owns JSON-RPC,
UTF-8 decoding and tool dispatch; no numerical dependency is preloaded here.
"""

from __future__ import annotations

from contextlib import contextmanager
import ctypes
from ctypes import wintypes
import io
import os
import sys
import threading


class _PipeReader(io.RawIOBase):
    """Borrow a pipe handle without closing it or changing its operating mode."""

    def __init__(self, handle, kernel, descriptor: int):
        super().__init__()
        self._handle = handle
        self._kernel = kernel
        self._descriptor = descriptor
        self._stop = threading.Event()

    def readable(self):
        return True

    def fileno(self):
        return self._descriptor

    def readinto(self, buffer):
        if self.closed:
            raise ValueError("read from closed MCP stdin")
        if not len(buffer):
            return 0
        while not self._stop.is_set():
            available = wintypes.DWORD()
            if not self._kernel.PeekNamedPipe(
                self._handle, None, 0, None, ctypes.byref(available), None
            ):
                return self._eof_or_raise()
            if not available.value:
                self._stop.wait(0.01)
                continue
            # This transport is the sole pipe reader. Never request missing data.
            count = min(len(buffer), available.value, io.DEFAULT_BUFFER_SIZE)
            target = (ctypes.c_char * count).from_buffer(buffer)
            received = wintypes.DWORD()
            if not self._kernel.ReadFile(
                self._handle, target, count, ctypes.byref(received), None
            ):
                return self._eof_or_raise()
            if received.value > count:
                raise OSError("MCP pipe returned an invalid byte count")
            return received.value
        return 0

    @staticmethod
    def _eof_or_raise():
        error = ctypes.get_last_error()
        if error in (109, 232, 233):  # Broken or disconnected pipe.
            return 0
        raise OSError(error, "MCP standard-input pipe read failed")

    def close(self):
        self._stop.set()
        super().close()  # The standard OS handle and CRT descriptor are borrowed.


def _reader_for_stdin(stream):
    if os.name != "nt":
        return None
    import msvcrt

    try:
        descriptor = stream.fileno()
    except (AttributeError, io.UnsupportedOperation):
        return None
    handle = msvcrt.get_osfhandle(descriptor)
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetFileType.argtypes = [wintypes.HANDLE]
    kernel.GetFileType.restype = wintypes.DWORD
    if kernel.GetFileType(handle) != 3:  # Only piped stdin needs this adapter.
        return None
    kernel.PeekNamedPipe.argtypes = [
        wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD, ctypes.c_void_p,
        ctypes.POINTER(wintypes.DWORD), ctypes.c_void_p,
    ]
    kernel.PeekNamedPipe.restype = wintypes.BOOL
    kernel.ReadFile.argtypes = [
        wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD,
        ctypes.POINTER(wintypes.DWORD), ctypes.c_void_p,
    ]
    kernel.ReadFile.restype = wintypes.BOOL
    return _PipeReader(handle, kernel, descriptor)


@contextmanager
def _pipe_stdin():
    previous = sys.stdin
    reader = _reader_for_stdin(previous)
    if reader is None:
        yield
        return
    replacement = io.TextIOWrapper(
        io.BufferedReader(reader), encoding="utf-8", errors="replace"
    )
    try:
        sys.stdin = replacement
        yield
    finally:
        sys.stdin = previous
        reader.close()  # Wake a pending read before closing the buffered wrapper.
        replacement.close()


def run_stdio(server):
    """Run the existing SDK transport with a scoped Windows pipe reader."""
    with _pipe_stdin():
        server.run(transport="stdio")
