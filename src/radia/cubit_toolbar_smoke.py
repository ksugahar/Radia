"""Cold-start smoke test for the persistent Radia Export Cubit toolbar.

Unlike the static toolbar-package tests, this command starts the real Cubit
GUI and asks Cubit's embedded PySide6 runtime what is actually visible.  It
runs two fresh starts by default to catch first-start timing failures and
loss of the persisted Coreform WorkflowToolbar on the next start.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SCHEMA = "radia.cubit-toolbar-smoke.v1"
PROBE_SCHEMA = "radia.cubit-toolbar-probe.v1"
EXPECTED_ACTIONS = [
    "Netgen Vol (.vol)",
    "GMSH (.msh)",
    "Nastran (.bdf)",
    "VTK (.vtk)",
    "FEMEEM",
    "MEG (ELF/MAGIC)",
]


def _find_cubit_exe() -> Path | None:
    explicit = os.environ.get("CUBIT_EXE")
    if explicit and Path(explicit).is_file():
        return Path(explicit)
    if sys.platform == "win32":
        program_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        candidates = sorted(
            glob.glob(str(Path(program_files) / "Coreform Cubit *" / "bin"
                          / "coreform_cubit.exe")),
            reverse=True,
        )
        if candidates:
            return Path(candidates[0])
    for name in ("coreform_cubit", "cubit"):
        found = shutil.which(name)
        if found:
            return Path(found)
    return None


def _cubit_pids() -> set[int]:
    """Return running Cubit process IDs without adding a psutil dependency."""
    if sys.platform != "win32":
        return set()
    result = subprocess.run(
        ["tasklist", "/FI", "IMAGENAME eq coreform_cubit.exe", "/FO", "CSV",
         "/NH"],
        capture_output=True,
        check=False,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    pids = set()
    for row in csv.reader(result.stdout.splitlines()):
        if len(row) < 2 or row[0].lower() != "coreform_cubit.exe":
            continue
        try:
            pids.add(int(row[1]))
        except ValueError:
            continue
    return pids


def _terminate_pids(pids: set[int]) -> None:
    for pid in sorted(pids):
        subprocess.run(
            ["taskkill", "/PID", str(pid), "/T", "/F"],
            capture_output=True,
            check=False,
        )


def validate_probe_result(payload: dict) -> list[str]:
    """Return display-contract violations from one Cubit GUI snapshot."""
    issues = []
    if payload.get("schema") != PROBE_SCHEMA:
        issues.append(f"unexpected probe schema: {payload.get('schema')!r}")
    if payload.get("main_window_visible") is not True:
        issues.append("Cubit main window is not visible")
    if payload.get("toolbar_count") != 1:
        issues.append(
            f"expected exactly one Radia Export toolbar, found "
            f"{payload.get('toolbar_count')!r}"
        )
    if payload.get("toolbar_visible") is not True:
        issues.append("Radia Export toolbar is not visible")
    if payload.get("toolbar_visible_region_nonempty") is not True:
        issues.append("Radia Export toolbar has no visible screen region")
    size = payload.get("toolbar_size")
    if not (
        isinstance(size, list) and len(size) == 2
        and all(isinstance(value, int) and value > 0 for value in size)
    ):
        issues.append(f"Radia Export toolbar has invalid size: {size!r}")

    actions = payload.get("toolbar_actions")
    if actions != EXPECTED_ACTIONS:
        issues.append(
            f"toolbar actions differ: expected {EXPECTED_ACTIONS!r}, got "
            f"{actions!r}"
        )
    visible = payload.get("action_visible")
    enabled = payload.get("action_enabled")
    if not isinstance(visible, dict):
        issues.append("per-action visibility report is missing")
    else:
        for name in EXPECTED_ACTIONS:
            if visible.get(name) is not True:
                issues.append(f"toolbar action is not visible: {name}")
    if not isinstance(enabled, dict):
        issues.append("per-action enabled report is missing")
    else:
        for name in EXPECTED_ACTIONS:
            if enabled.get(name) is not True:
                issues.append(f"toolbar action is not enabled: {name}")
    if payload.get("toolbar_menu_has_radia_export") is not True:
        issues.append("Radia Export is absent from Cubit's toolbar menu")
    if payload.get("unsupported_top_level_menu_present") is not False:
        issues.append("unsupported top-level Radia Export QMenu is present")
    if payload.get("ok") is not True:
        issues.append("Cubit-side display probe did not report success")
    return list(dict.fromkeys(issues))


def _probe_path() -> Path:
    return Path(__file__).resolve().parent / "panels" / "cubit_toolbar_probe.py"


def _write_bootstrap(path: Path, probe_path: Path) -> None:
    probe = str(probe_path.resolve())
    source = (
        "#!python\n"
        f"exec(compile(open({probe!r}, encoding='utf-8').read(), "
        f"{probe!r}, 'exec'))\n"
    )
    path.write_text(source, encoding="utf-8")


def _run_one(cubit_exe: Path, work: Path, timeout: float) -> dict:
    result_path = work / "result.json"
    bootstrap_path = work / "probe_bootstrap.py"
    probe_path = _probe_path()
    if not probe_path.is_file():
        return {"ok": False, "issues": [f"probe script missing: {probe_path}"]}
    _write_bootstrap(bootstrap_path, probe_path)

    before = _cubit_pids()
    if before:
        return {
            "ok": False,
            "issues": [
                (
                    "Cubit is already running; close it before the cold-start "
                    f"test (pids={sorted(before)})"
                )
            ],
        }

    env = os.environ.copy()
    env["RADIA_TOOLBAR_PROBE_RESULT"] = str(result_path)
    env["RADIA_TOOLBAR_PROBE_TIMEOUT"] = str(max(5.0, timeout - 10.0))
    launcher = subprocess.Popen(
        [str(cubit_exe), "-nojournal", str(bootstrap_path)],
        cwd=str(work),
        env=env,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )

    deadline = time.monotonic() + timeout
    spawned = set()
    while time.monotonic() < deadline and not result_path.is_file():
        spawned.update(_cubit_pids() - before)
        time.sleep(0.1)

    if not result_path.is_file():
        spawned.update(_cubit_pids() - before)
        _terminate_pids(spawned)
        return {
            "ok": False,
            "issues": [
                f"Cubit GUI produced no toolbar result within {timeout:.1f}s"
            ],
            "launcher_returncode": launcher.poll(),
        }

    try:
        payload = json.loads(result_path.read_text(encoding="utf-8-sig"))
    except (OSError, json.JSONDecodeError) as exc:
        spawned.update(_cubit_pids() - before)
        _terminate_pids(spawned)
        return {"ok": False, "issues": [f"invalid probe result: {exc}"]}

    issues = validate_probe_result(payload)
    close_deadline = time.monotonic() + 15.0
    while time.monotonic() < close_deadline:
        spawned.update(_cubit_pids() - before)
        alive = spawned & _cubit_pids()
        if not alive:
            break
        time.sleep(0.1)
    else:
        alive = spawned & _cubit_pids()
        _terminate_pids(alive)
        issues.append(f"Cubit GUI did not close after probe (pids={sorted(alive)})")

    return {
        "ok": not issues,
        "issues": issues,
        "probe": payload,
        "launcher_returncode": launcher.poll(),
    }


def run_smoke_test(*, restarts: int = 2, timeout: float = 45.0,
                   keep: bool = False, report_json: str = "") -> int:
    if restarts < 1:
        print("[FAIL] --restarts must be at least 1")
        return 1
    cubit_exe = _find_cubit_exe()
    if cubit_exe is None:
        print("[FAIL] Coreform Cubit executable not found; set CUBIT_EXE")
        return 1
    if _cubit_pids():
        print("[FAIL] Cubit is already running. Close it before this cold-start test.")
        return 1

    temp_root = Path(r"C:\temp") if os.name == "nt" else Path(
        tempfile.gettempdir()
    )
    temp_root.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="cubit-toolbar-smoke-", dir=temp_root))
    print(f"  Cubit: {cubit_exe}")
    print(f"  Work:  {root}")

    runs = []
    for index in range(1, restarts + 1):
        work = root / f"start-{index}"
        work.mkdir()
        print(f"  Cold start {index}/{restarts} ...", end="", flush=True)
        result = _run_one(cubit_exe, work, timeout)
        runs.append(result)
        if result.get("ok"):
            probe = result["probe"]
            print(
                f" OK ({probe.get('elapsed_seconds')}s, "
                f"{len(probe.get('toolbar_actions', []))} actions)"
            )
        else:
            print(" FAIL")
            for issue in result.get("issues", []):
                print(f"       - {issue}")
            break

    report = {
        "schema": SCHEMA,
        "passed": len(runs) == restarts and all(run.get("ok") for run in runs),
        "cubit_executable": str(cubit_exe),
        "requested_restarts": restarts,
        "completed_restarts": len(runs),
        "runs": runs,
    }
    destination = Path(report_json) if report_json else root / "report.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )

    if report["passed"]:
        print(
            "[OK] Radia Export is visible and complete on every cold start: "
            "one toolbar, six visible/enabled actions, persistent toolbar-menu "
            "entry, no unsupported top-level QMenu."
        )
        if not keep and not report_json:
            shutil.rmtree(root, ignore_errors=True)
        else:
            print(f"  Report: {destination}")
        return 0

    print(f"[FAIL] Radia Export GUI display contract failed. Report: {destination}")
    print(f"       workdir retained for diagnosis: {root}")
    return 1


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="cubit-toolbar-smoke-test",
        description=(
            "Cold-start real Coreform Cubit and verify that the persistent "
            "Radia Export toolbar is actually visible and complete."
        ),
    )
    parser.add_argument(
        "--restarts", type=int, default=2,
        help="number of independent GUI cold starts (default: 2)",
    )
    parser.add_argument(
        "--timeout", type=float, default=45.0,
        help="seconds allowed for each GUI start (default: 45)",
    )
    parser.add_argument("--keep", action="store_true")
    parser.add_argument(
        "--report-json", default="",
        help="optional persistent JSON report path",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("  Radia Export: real Cubit GUI cold-start smoke test")
    print("=" * 60)
    raise SystemExit(run_smoke_test(
        restarts=args.restarts,
        timeout=args.timeout,
        keep=args.keep,
        report_json=args.report_json,
    ))


if __name__ == "__main__":
    main()
