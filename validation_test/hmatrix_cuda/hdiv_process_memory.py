"""Windows process samples; not global memory or a guaranteed phase peak."""
import ctypes
import threading


class ProcessMemory:
    def __init__(self,pid):
        self.pid=pid;self.samples=[];self.errors=[];self.stop=threading.Event()

    def __enter__(self):
        size=ctypes.c_size_t
        class Counters(ctypes.Structure):
            _fields_=[('cb',ctypes.c_ulong),('faults',ctypes.c_ulong)]+[(name,size) for name in ('peak_rss','rss','peak_paged','paged','peak_nonpaged','nonpaged','pagefile','peak_pagefile','private')]
        self.counters= Counters
        self.kernel=ctypes.WinDLL('kernel32',use_last_error=True)
        self.kernel.OpenProcess.argtypes=[ctypes.c_ulong,ctypes.c_int,ctypes.c_ulong]
        self.kernel.OpenProcess.restype=ctypes.c_void_p
        self.kernel.CloseHandle.argtypes=[ctypes.c_void_p]
        self.psapi=ctypes.WinDLL('psapi',use_last_error=True)
        self.psapi.GetProcessMemoryInfo.argtypes=[ctypes.c_void_p,ctypes.POINTER(Counters),ctypes.c_ulong]
        self.handle=self.kernel.OpenProcess(0x410,False,self.pid)
        if not self.handle:raise ctypes.WinError(ctypes.get_last_error())
        try:self.sample()
        except BaseException:
            self.kernel.CloseHandle(self.handle);raise
        def observe():
            try:
                while not self.stop.wait(.05):self.sample()
            except BaseException as exc:self.errors.append(repr(exc))
        self.thread=threading.Thread(target=observe,daemon=True)
        self.thread.start()
        return self

    def sample(self):
        values=self.counters();values.cb=ctypes.sizeof(values)
        if not self.psapi.GetProcessMemoryInfo(self.handle,ctypes.byref(values),values.cb):
            raise ctypes.WinError(ctypes.get_last_error())
        self.samples.append(dict(rss_bytes=values.rss,private_bytes=values.private,lifetime_peak_rss_bytes=values.peak_rss))

    def __exit__(self,*args):
        self.stop.set();self.thread.join()
        try:self.sample()
        finally:self.kernel.CloseHandle(self.handle)
        if self.errors:raise RuntimeError('Memory sampling failed: '+str(self.errors))

    def result(self):
        return dict(before=self.samples[0],after=self.samples[-1],sample_count=len(self.samples),period_seconds=.05,
                    sampled_peak_rss_bytes=max(s['rss_bytes'] for s in self.samples),
                    sampled_peak_private_bytes=max(s['private_bytes'] for s in self.samples),
                    limitation='Process samples include baseline and imported runtimes, may miss short spikes; lifetime RSS peak is not a phase peak; no global host/device peak')
