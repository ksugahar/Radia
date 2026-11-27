"""Cold-start smoke test for the persistent Netgen Export Cubit toolbar.

Unlike the static toolbar-package tests, this command starts the real Cubit
GUI and asks Cubit's embedded PySide6 runtime what is actually visible.  It
runs two fresh starts by default to catch first-start timing failures and
loss of the persisted Coreform WorkflowToolbar on the next start.

This is the human-facing GUI release test of the exporter toolbar; LLM/MCP
Cubit execution never runs it.  Run it as
``python -m radia_mcp.cubit.toolbar_smoke`` after
``python -m radia_mcp.cubit.install``.
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

from radia_mcp.cubit.layout import is_supported_cubit_root, runtime_environment

SCHEMA = "radia-mcp.cubit-toolbar-smoke.v1"
PROBE_SCHEMA = "radia-mcp.cubit-toolbar-probe.v1"
TOOLBAR_TITLE = "Netgen Export"
EXPECTED_ACTIONS = [
    "Export Netgen (.vol)...",
]


def _find_cubit_exe() -> Path | None:
    explicit = os.environ.get("CUBIT_EXE")
    if explicit:
        path = Path(explicit)
        return path if path.is_file() and is_supported_cubit_root(path.resolve().parent.parent) else None
    if sys.platform == "win32":
        program_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        candidates = sorted(
            glob.glob(str(Path(program_files) / "Coreform Cubit *" / "bin"
                          / "clarox.exe")) +
            glob.glob(str(Path(program_files) / "Coreform Cubit *" / "bin"
                          / "coreform_cubit.exe")),
            reverse=True,
        )
        for candidate in candidates:
            path = Path(candidate)
            if is_supported_cubit_root(path.parent.parent):
                return path
    for name in ("coreform_cubit", "cubit"):
        found = shutil.which(name)
        if found and is_supported_cubit_root(Path(found).resolve().parent.parent):
            return Path(found)
    return None


def _cubit_pids() -> set[int]:
    """Return running Cubit process IDs without adding a psutil dependency."""
    if sys.platform != "win32":
        return set()
    result = subprocess.run(
        ["tasklist", "/FO", "CSV", "/NH"],
        capture_output=True,
        check=False,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    pids = set()
    for row in csv.reader(result.stdout.splitlines()):
        if len(row) < 2 or row[0].lower() not in {"coreform_cubit.exe", "clarox.exe", "cubitx.exe"}:
            continue
        try:
            pids.add(int(row[1]))
        except ValueError:
            continue
    return pids


def _terminate_owned(launcher) -> None:
    """Never infer ownership from newly appearing, potentially user-owned PIDs."""
    if launcher.poll() is None:
        launcher.terminate()
        try:
            launcher.wait(timeout=10)
        except subprocess.TimeoutExpired:
            launcher.kill()
            launcher.wait(timeout=10)


def validate_probe_result(payload: dict) -> list[str]:
    """Return display-contract violations from one Cubit GUI snapshot."""
    issues = []
    if payload.get("schema") != PROBE_SCHEMA:
        issues.append(f"unexpected probe schema: {payload.get('schema')!r}")
    if payload.get("main_window_visible") is not True:
        issues.append("Cubit main window is not visible")
    if payload.get("toolbar_count") != 1:
        issues.append(
            f"expected exactly one {TOOLBAR_TITLE} toolbar, found "
            f"{payload.get('toolbar_count')!r}"
        )
    if payload.get("toolbar_visible") is not True:
        issues.append(f"{TOOLBAR_TITLE} toolbar is not visible")
    if payload.get("toolbar_visible_region_nonempty") is not True:
        issues.append(f"{TOOLBAR_TITLE} toolbar has no visible screen region")
    size = payload.get("toolbar_size")
    if not (
        isinstance(size, list) and len(size) == 2
        and all(isinstance(value, int) and value > 0 for value in size)
    ):
        issues.append(f"{TOOLBAR_TITLE} toolbar has invalid size: {size!r}")

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
    if payload.get("toolbar_menu_has_entry") is not True:
        issues.append(f"{TOOLBAR_TITLE} is absent from Cubit's toolbar menu")
    if payload.get("toolbar_owner") != "WorkflowToolbar":
        issues.append("toolbar is not owned by Coreform WorkflowToolbar")
    if payload.get("ok") is not True:
        issues.append("Cubit-side display probe did not report success")
    return list(dict.fromkeys(issues))


def _probe_path() -> Path:
    return Path(__file__).resolve().parent / "gui" / "toolbar_probe.py"


def _write_bootstrap(path: Path, probe_path: Path) -> None:
    probe = str(probe_path.resolve())
    source = (
        "#!python\n"
        f"exec(compile(open({probe!r}, encoding='utf-8').read(), "
        f"{probe!r}, 'exec'), {{'__name__': '__main__'}})\n"
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

    env = runtime_environment(cubit_exe.parent)
    env["RADIA_MCP_CUBIT_TOOLBAR_PROBE_RESULT"] = str(result_path)
    env["RADIA_MCP_CUBIT_TOOLBAR_PROBE_TIMEOUT"] = str(max(5.0, timeout - 10.0))
    modern = cubit_exe.name.lower() == "clarox.exe"
    command = [str(cubit_exe)]
    if modern:
        command += ["-prompt", "Cubit"]
    command += ["-nojournal"]
    if not modern:
        command += ["-commandplugindir", str(cubit_exe.parent / "plugins")]
    command.append(str(bootstrap_path))
    with (work / "launcher.log").open("wb") as output:
        launcher = subprocess.Popen(
            command,
            cwd=str(work),
            env=env,
            stdin=subprocess.PIPE,
            stdout=output,
            stderr=subprocess.STDOUT,
        )

    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and not result_path.is_file():
        if launcher.poll() is not None:
            break
        time.sleep(0.1)

    if not result_path.is_file():
        _terminate_owned(launcher)
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
        _terminate_owned(launcher)
        return {"ok": False, "issues": [f"invalid probe result: {exc}"]}

    issues = validate_probe_result(payload)
    try:
        returncode = launcher.wait(timeout=15)
    except subprocess.TimeoutExpired:
        _terminate_owned(launcher)
        issues.append("Cubit GUI did not close after probe")
    else:
        if returncode != 0:
            issues.append(f"Cubit GUI exited with failure: {returncode}")

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
    root = Path(tempfile.mkdtemp(prefix="radia-mcp-toolbar-smoke-", dir=temp_root))
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
            f"[OK] {TOOLBAR_TITLE} is visible and complete on every cold start: "
            f"one toolbar, {len(EXPECTED_ACTIONS)} visible/enabled action(s), "
            "persistent toolbar-menu entry, Coreform-owned WorkflowToolbar."
        )
        if not keep and not report_json:
            shutil.rmtree(root, ignore_errors=True)
        else:
            print(f"  Report: {destination}")
        return 0

    print(f"[FAIL] {TOOLBAR_TITLE} GUI display contract failed. Report: {destination}")
    print(f"       workdir retained for diagnosis: {root}")
    return 1


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="python -m radia_mcp.cubit.toolbar_smoke",
        description=(
            "Cold-start real Coreform Cubit and verify that the persistent "
            f"{TOOLBAR_TITLE} toolbar is actually visible and complete."
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
    args = parser.parse_args(argv)

    print("=" * 60)
    print(f"  {TOOLBAR_TITLE}: real Cubit GUI cold-start smoke test")
    print("=" * 60)
    raise SystemExit(run_smoke_test(
        restarts=args.restarts,
        timeout=args.timeout,
        keep=args.keep,
        report_json=args.report_json,
    ))


if __name__ == "__main__":
    main()
