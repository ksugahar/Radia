"""Windows-owned worker trees; GUI lifetime is deliberately outside these jobs."""
from __future__ import annotations

import asyncio
import ctypes as c
from ctypes import wintypes as w
import os
import subprocess
import sys
import time


class _Job:
    def __init__(self):
        if sys.platform != 'win32':
            raise RuntimeError('Mesh viewer worker ownership requires Windows')
        self.k = k = c.WinDLL('kernel32', use_last_error=True)
        class Basic(c.Structure):
            _fields_ = [('user', c.c_longlong), ('kernel', c.c_longlong), ('flags', w.DWORD),
                        ('min_ws', c.c_size_t), ('max_ws', c.c_size_t), ('active', w.DWORD),
                        ('affinity', c.c_size_t), ('priority', w.DWORD), ('scheduling', w.DWORD)]
        class IO(c.Structure):
            _fields_ = [(name, c.c_ulonglong) for name in ('ro', 'wo', 'oo', 'rb', 'wb', 'ob')]
        class Limits(c.Structure):
            _fields_ = [('basic', Basic), ('io', IO), ('pm', c.c_size_t), ('jm', c.c_size_t),
                        ('peak_p', c.c_size_t), ('peak_j', c.c_size_t)]
        k.CreateJobObjectW.argtypes = [w.LPVOID, w.LPCWSTR]; k.CreateJobObjectW.restype = w.HANDLE
        k.SetInformationJobObject.argtypes = [w.HANDLE, c.c_int, w.LPVOID, w.DWORD]
        k.QueryInformationJobObject.argtypes = [w.HANDLE, c.c_int, w.LPVOID, w.DWORD, w.LPVOID]
        k.OpenProcess.argtypes = [w.DWORD, w.BOOL, w.DWORD]; k.OpenProcess.restype = w.HANDLE
        k.AssignProcessToJobObject.argtypes = [w.HANDLE, w.HANDLE]
        k.IsProcessInJob.argtypes = [w.HANDLE, w.HANDLE, c.POINTER(w.BOOL)]
        k.CloseHandle.argtypes = [w.HANDLE]
        k.TerminateJobObject.argtypes = [w.HANDLE, w.UINT]
        self.handle = k.CreateJobObjectW(None, None)
        self.check(self.handle)
        try:
            limits = Limits(); limits.basic.flags = 0x2000  # KILL_ON_JOB_CLOSE, no breakaway.
            self.check(k.SetInformationJobObject(self.handle, 9, c.byref(limits), c.sizeof(limits)))
        except BaseException:
            k.CloseHandle(self.handle)
            raise

    @staticmethod
    def check(value):
        if not value: raise c.WinError(c.get_last_error())

    def attach_resume(self, pid):
        k = self.k
        process = k.OpenProcess(0x100000 | 0x1000 | 0x100 | 1, False, pid)
        self.check(process)
        try:
            self.check(k.AssignProcessToJobObject(self.handle, process))
            inside = w.BOOL()
            self.check(k.IsProcessInJob(process, self.handle, c.byref(inside)))
            if not inside.value: raise RuntimeError('Worker not in its explicit job')
            class Thread(c.Structure):
                _fields_ = [(name, w.DWORD) for name in ('size', 'usage', 'tid', 'pid')] + [
                    ('base', w.LONG), ('delta', w.LONG), ('flags', w.DWORD)]
            k.CreateToolhelp32Snapshot.argtypes = [w.DWORD, w.DWORD]
            k.CreateToolhelp32Snapshot.restype = w.HANDLE
            k.Thread32First.argtypes = [w.HANDLE, c.POINTER(Thread)]
            k.Thread32Next.argtypes = [w.HANDLE, c.POINTER(Thread)]
            k.OpenThread.argtypes = [w.DWORD, w.BOOL, w.DWORD]; k.OpenThread.restype = w.HANDLE
            k.ResumeThread.argtypes = [w.HANDLE]; k.ResumeThread.restype = w.DWORD
            snapshot = k.CreateToolhelp32Snapshot(4, 0)
            if snapshot == c.c_void_p(-1).value: raise c.WinError(c.get_last_error())
            try:
                row = Thread(); row.size = c.sizeof(row); resumed = 0
                more = k.Thread32First(snapshot, c.byref(row))
                while more:
                    if row.pid == pid:
                        thread = k.OpenThread(2, False, row.tid); self.check(thread)
                        try:
                            if k.ResumeThread(thread) == 0xffffffff: raise c.WinError(c.get_last_error())
                            resumed += 1
                        finally: k.CloseHandle(thread)
                    more = k.Thread32Next(snapshot, c.byref(row))
                if not resumed: raise RuntimeError('Suspended worker has no owned main thread')
            finally: k.CloseHandle(snapshot)
        finally: k.CloseHandle(process)

    def pids(self):
        class Pids(c.Structure):
            _fields_ = [('assigned', w.DWORD), ('count', w.DWORD), ('ids', c.c_size_t * 256)]
        row = Pids()
        self.check(self.k.QueryInformationJobObject(self.handle, 3, c.byref(row), c.sizeof(row), None))
        if row.count > 256: raise RuntimeError('Worker tree exceeds ownership-query capacity')
        return list(row.ids[:row.count])

    def terminate(self): self.check(self.k.TerminateJobObject(self.handle, 1))
    def close(self): self.k.CloseHandle(self.handle)


async def _drain(stream):
    data = bytearray()
    while chunk := await stream.read(65536):
        data.extend(chunk)
        if len(data) > 1024 * 1024:
            raise RuntimeError('Worker diagnostics exceed the 1 MiB capture limit')
    return bytes(data)


async def run_worker(argv, *, env, deadline, phase):
    """Await one worker plus descendants within an absolute monotonic deadline."""
    remaining = deadline - time.monotonic()
    if remaining <= 5:
        raise RuntimeError(f'{phase}: transaction deadline exhausted before worker launch')
    job, proc, readers = _Job(), None, []
    spawn = None
    try:
        child_env = {**env, 'PYTHONIOENCODING': 'utf-8', 'PYTHONUTF8': '1',
                     'PYTHONDONTWRITEBYTECODE': '1'}
        spawn = asyncio.create_task(asyncio.create_subprocess_exec(
            *argv, env=child_env, stdin=subprocess.DEVNULL, stdout=asyncio.subprocess.PIPE,
            stderr=asyncio.subprocess.PIPE,
            creationflags=0x4 | subprocess.CREATE_NO_WINDOW))  # CREATE_SUSPENDED
        proc = await asyncio.shield(spawn)
        job.attach_resume(proc.pid)
        readers = [asyncio.create_task(_drain(proc.stdout)), asyncio.create_task(_drain(proc.stderr))]
        async def complete():
            while proc.returncode is None:
                for reader in readers:
                    if reader.done() and reader.exception(): raise reader.exception()
                await asyncio.sleep(.01)
            settle = time.monotonic() + .2
            while job.pids() and time.monotonic() < settle:
                await asyncio.sleep(.01)
            if job.pids():
                raise RuntimeError(f'{phase}: worker exited with live descendants')
            stdout, stderr = await asyncio.gather(*readers)
            return subprocess.CompletedProcess(argv, proc.returncode, stdout, stderr)
        try:
            result = await asyncio.wait_for(complete(), deadline - time.monotonic() - 5)
        except TimeoutError as error:
            raise RuntimeError(f'{phase}: transaction timed out (remaining budget {remaining:.3f}s); owned worker tree terminated') from error
        if result.returncode:
            raise RuntimeError(f'{phase}: ' + result.stderr.decode('utf-8', errors='replace').strip())
        return result
    finally:
        try:
            if proc is None and spawn is not None:
                try: proc = await asyncio.wait_for(asyncio.shield(spawn), 3)
                except TimeoutError:
                    spawn.cancel()
                    await asyncio.gather(spawn, return_exceptions=True)
            job.terminate()
            if proc is not None and proc.returncode is None:
                proc.kill()  # Also covers assignment failure: this exact suspended launcher only.
            for reader in readers:
                if not reader.done(): reader.cancel()
            if readers: await asyncio.gather(*readers, return_exceptions=True)
            if proc is not None:
                async def discard(stream):
                    while await stream.read(65536): pass
                await asyncio.wait_for(asyncio.gather(proc.wait(), discard(proc.stdout),
                                                      discard(proc.stderr)), 3)
            end = time.monotonic() + 2
            while job.pids():
                if time.monotonic() > end: raise RuntimeError(f'{phase}: owned tree failed to terminate')
                await asyncio.sleep(.01)
        finally:
            job.close()
            for reader in readers:
                if not reader.done(): reader.cancel()
            if readers: await asyncio.gather(*readers, return_exceptions=True)


def transaction_deadline():
    import math
    try: budget = float(os.environ.get('RADIA_MESH_VIEWER_TIMEOUT_SECONDS', '300'))
    except ValueError as error: raise ValueError('Viewer timeout must be a finite 30..900 seconds') from error
    if not math.isfinite(budget) or not 30 <= budget <= 900:
        raise ValueError('Viewer timeout must be a finite 30..900 seconds')
    return time.monotonic() + budget
