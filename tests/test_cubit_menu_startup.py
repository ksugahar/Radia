"""Run the menu-persistence regression in Cubit's own PySide6 runtime."""

from __future__ import annotations

import glob
import json
import os
import subprocess
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "tests" / "cubit_menu_runtime_probe.py"
EXPORT_MENU = ROOT / "src" / "radia" / "panels" / "radia_export_menu.py"


def _find_cubit_python() -> Path | None:
    explicit = os.environ.get("CUBIT_PYTHON")
    if explicit and Path(explicit).is_file():
        return Path(explicit)

    cubit_path = os.environ.get("CUBIT_PATH")
    if cubit_path:
        root = Path(cubit_path)
        bin_dir = root if root.name.lower() == "bin" else root / "bin"
        candidate = bin_dir / "python3" / "python.exe"
        if candidate.is_file():
            return candidate

    if os.name == "nt":
        program_files = os.environ.get("ProgramFiles", r"C:\Program Files")
        candidates = sorted(
            glob.glob(str(Path(program_files) / "Coreform Cubit *" / "bin"
                           / "python3" / "python.exe")),
            reverse=True,
        )
        if candidates:
            return Path(candidates[0])
    return None


def test_cubit_runtime_restores_menu_after_repeated_stock_menu_rebuilds():
    cubit_python = _find_cubit_python()
    if cubit_python is None:
        pytest.skip("Coreform Cubit embedded Python is not installed")

    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    proc = subprocess.run(
        [str(cubit_python), str(PROBE), str(EXPORT_MENU)],
        cwd=str(ROOT),
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    json_line = next(
        (line for line in reversed(proc.stdout.splitlines())
         if line.lstrip().startswith("{")),
        "",
    )
    assert proc.returncode == 0, (
        f"Cubit PySide6 menu probe failed (rc={proc.returncode})\n"
        f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
    )
    assert json_line, f"probe returned no JSON payload:\n{proc.stdout}"
    result = json.loads(json_line)

    assert result["ok"] is True
    assert result["runtime"]["python"].startswith("3.10")
    assert result["observer"]["installed_before_main_window"] is True
    assert result["observer"]["survived_gc"] is True
    assert result["observer"]["parented_to_qapplication"] is True
    assert result["replay_menu_count"] == 1
    assert len(result["rebuilds"]) == 4
    assert all(row["menu_count"] == 1 for row in result["rebuilds"])
    assert all(row["action_count"] == 6 for row in result["rebuilds"])
