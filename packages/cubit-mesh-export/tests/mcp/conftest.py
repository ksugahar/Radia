"""Fail on surviving test-owned Cubit runners; never terminate user processes."""
from pathlib import PurePosixPath
import os
import subprocess
import time

import pytest


def _is_cubit_runner(argv):
    """Match script/module arguments, not arbitrary text mentioning a runner."""
    args = [str(arg).replace("\\", "/") for arg in argv or []]
    modules = {"cubit_mesh_export.mcp.daemon", "cubit_mesh_export.mcp.bootstrap"}
    if any(arg == "-m" and index + 1 < len(args) and args[index + 1] in modules
           for index, arg in enumerate(args)):
        return True
    return any(PurePosixPath(arg).parts[-3:] in {
        ("cubit_mesh_export", "mcp", "daemon.py"),
        ("cubit_mesh_export", "mcp", "bootstrap.py"),
    } for arg in args[1:])


def _process_snapshot(psutil):
    rows = {}
    for proc in psutil.process_iter(["pid", "ppid", "create_time", "cmdline"]):
        try:
            row = proc.info
            if row["create_time"] is not None:
                rows[row["pid"]] = dict(row)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return rows


def _new_owned_runners(before, after, owners):
    """Use process identity and ancestry, including an exited recorded parent."""
    leaked = {}
    for pid, row in after.items():
        identity = (pid, row["create_time"])
        if identity in before or not _is_cubit_runner(row["cmdline"]):
            continue
        current = row
        seen = set()
        while current["pid"] not in seen:
            current_pid = current["pid"]
            seen.add(current_pid)
            if owners.get(current_pid) == current["create_time"]:
                leaked[pid] = row
                break
            parent = current["ppid"]
            if parent in after:
                current = after[parent]
            elif parent in owners and owners[parent] <= current["create_time"]:
                # A test-spawned MCP parent exited, leaving its daemon alive.
                leaked[pid] = row
                break
            else:
                break
    return leaked


@pytest.fixture(scope="session", autouse=True)
def cubit_process_leak_gate():
    try:
        import psutil
    except ImportError:
        pytest.fail("Cubit MCP lifetime tests require psutil; install the test dependency.")
    initial = _process_snapshot(psutil)
    before = {(pid, row["create_time"]) for pid, row in initial.items()}
    owners = {os.getpid(): psutil.Process().create_time()}
    original_init = subprocess.Popen.__init__

    def record_child(proc, *args, **kwargs):
        started = time.time()
        original_init(proc, *args, **kwargs)
        try:
            owners[proc.pid] = psutil.Process(proc.pid).create_time()
        except psutil.NoSuchProcess:
            # Retain the parent PID even if it exits before inspection.
            owners[proc.pid] = started

    # Patch the initializer, not the class: preserve isinstance/context-manager
    # behavior and allow individual tests to install their normal Popen fakes.
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(subprocess.Popen, "__init__", record_child)
        yield
    leaked = _new_owned_runners(before, _process_snapshot(psutil), dict(owners))
    if leaked:
        detail = "\n".join(f"  pid {pid}: {' '.join(row['cmdline'])[:160]}"
                           for pid, row in leaked.items())
        pytest.fail("Test-owned Cubit runners survived teardown:\n" + detail +
                    "\nClose owned CubitSession instances in test cleanup. "
                    "No process was terminated by this gate.", pytrace=False)
