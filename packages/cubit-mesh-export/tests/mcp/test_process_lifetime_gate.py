"""Synthetic process tables test ownership without launching or killing Cubit."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys

import pytest

spec = importlib.util.spec_from_file_location(
    "cubit_lifetime_gate", Path(__file__).with_name("conftest.py"))
gate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(gate)


@pytest.mark.parametrize("argv, expected", [
    (["python", "C:/pkg/cubit_mesh_export/mcp/daemon.py"], True),
    (["python", r"C:\pkg\cubit_mesh_export\mcp\bootstrap.py"], True),
    (["python", "-m", "cubit_mesh_export.mcp.daemon"], True),
    (["coreform_cubit.exe", "-nojournal"], False),
    (["python", "C:/unrelated/daemon.py"], False),
    (["python", "-c", "print('cubit_mesh_export.mcp.daemon')"], False),
])
def test_runner_signature(argv, expected):
    assert gate._is_cubit_runner(argv) is expected


def row(pid, parent, created=10, runner=True):
    return dict(pid=pid, ppid=parent, create_time=created,
                cmdline=["python", "-m", "cubit_mesh_export.mcp.daemon"]
                if runner else ["python", "-m", "pytest"])


def test_only_new_test_owned_runners_fail():
    after = {1: row(1, 0, 1, False), 2: row(2, 1), 3: row(3, 1),
             4: row(4, 999), 5: row(5, 1, runner=False)}
    assert set(gate._new_owned_runners({(3, 10)}, after, {1: 1})) == {2}


def test_exited_recorded_parent_does_not_hide_orphan():
    assert set(gate._new_owned_runners(set(), {3: row(3, 2)}, {2: 5})) == {3}


def test_reused_parent_pid_and_unrelated_process_are_not_owned():
    after = {2: row(2, 99, 20, False), 3: row(3, 2, 21)}
    assert gate._new_owned_runners(set(), after, {2: 5}) == {}


def test_exited_runners_do_not_fail():
    assert gate._new_owned_runners({(2, 10)}, {}, {1: 1}) == {}


def test_real_python_child_is_detected_and_reaped_without_cubit(tmp_path):
    import psutil

    script = tmp_path / "cubit_mesh_export" / "mcp" / "daemon.py"
    script.parent.mkdir(parents=True)
    script.write_text("import sys\nsys.stdin.read()\n", encoding="utf-8")
    owners = {os.getpid(): psutil.Process().create_time()}
    with subprocess.Popen([sys.executable, str(script)], stdin=subprocess.PIPE) as proc:
        try:
            assert proc.pid in gate._new_owned_runners(
                set(), gate._process_snapshot(psutil), owners)
        finally:
            proc.communicate(timeout=5)  # EOF, not process termination.
    assert proc.pid not in gate._new_owned_runners(
        set(), gate._process_snapshot(psutil), owners)
