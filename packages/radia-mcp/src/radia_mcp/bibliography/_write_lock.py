"""Cooperative, fail-fast OS locks; nonparticipating editors need snapshot checks."""
from contextlib import contextmanager
import os
from pathlib import Path


@contextmanager
def target_lock(path: Path):
    # Keep the sidecar inode: unlinking it after unlock races a waiting opener.
    target = path.resolve()
    sidecar = target.with_name("." + target.name + ".radia-write-lock")
    with sidecar.open("a+b") as handle:
        if sidecar.stat().st_size == 0:
            handle.write(b"0")
            handle.flush()
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise OSError(f"bibliography target busy: {target}") from exc
        try:
            yield
        finally:
            handle.seek(0)
            if os.name == "nt":
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
