"""Sample process RSS at engine boundaries, never call it isolated engine memory."""
import functools
import threading
import time


def instrument(function, emit, read_rss, interval=0.05):
    if interval <= 0:
        raise ValueError('interval must be positive')

    @functools.wraps(function)
    def wrapped(*args, **kwargs):
        stopped = threading.Event()
        samples = []
        failures = []

        def sample():
            try:
                samples.append((time.perf_counter(), int(read_rss())))
            except Exception as exc:
                failures.append(repr(exc))

        def poll():
            while not stopped.wait(interval):
                sample()

        sample()
        worker = threading.Thread(target=poll, name='phase-rss', daemon=False)
        worker.start()
        completed = False
        try:
            result = function(*args, **kwargs)
            completed = True
            return result
        finally:
            stopped.set()
            worker.join()
            sample()
            emit('engine_phase_memory', engine=function.__name__, completed=completed,
                 rss_start_bytes=samples[0][1] if samples else None,
                 rss_end_bytes=samples[-1][1] if samples else None,
                 sampled_max_process_rss_bytes=max((x[1] for x in samples), default=None),
                 samples=len(samples), sampling_interval_s=interval, errors=failures,
                 memory_scope='Sampled whole-process RSS during this engine; may retain previous-engine allocations; sampling can miss peaks. Not isolated engine peak memory.')
    return wrapped
