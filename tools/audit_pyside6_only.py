#!/usr/bin/env python3
"""Audit: Cubit toolbar has no legacy Qt and normal Radia stays no-PySide.

Enforces the Radia policy "no Qt5 / PyQt5 fallbacks" while keeping the Cubit
mesh-export toolbar as the only PySide surface.  Normal Radia Python must not
depend on PySide6; the toolbar runs inside Coreform Cubit's embedded Python.
When PySide6 is installed in the current interpreter, we smoke-test the toolbar
headlessly.  When it is not installed, the smoke is skipped because the target
runtime is Cubit's embedded Python, not this interpreter.

Static checks:
  1. ZERO real PyQt5 / PySide2 / PyQt6 import statements in tracked *.py
     (lint-rule regexes and "no PyQt5 fallback" comments do NOT count --
     only actual ``import``/``from ... import`` lines).
  2. Cubit toolbar modules import PySide6.
  3. radia pyproject.toml does not declare a PySide6 dependency or retired
     standalone panel scripts.
  4. cubit_mesh_export.ccm (if present in the package source) has NO Qt5 DLL
     dependency (pefile import-table scan).

Headless check (isolated subprocess, offscreen Qt):
  5. If PySide6 is installed, ExportDialog builds all 6 formats emitting valid
     ``export ... overwrite`` commands.  If PySide6 is absent, this is a skip,
     not a failure.

Run:  python tools/audit_pyside6_only.py
Exit 0 = no legacy Qt and Cubit toolbar smoke healthy/skipped; non-zero =
issues (listed).

The machine/Cubit half (cubit-plugin-install --verify-only, cubit-smoke-test,
cross-machine SSH to 100号機 / mdx / hibino) lives in the release/deploy gates.
"""
from __future__ import annotations

import json
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

CUBIT_TOOLBAR_MODULES = [
    "src/radia/panels/radia_export_menu.py",
    "src/radia/panels/register_toolbar.py",
]


# ----------------------------------------------------------------------
# Static checks
# ----------------------------------------------------------------------
def check_no_legacy_qt_imports() -> list[str]:
    """git grep for real PyQt5/PySide2/PyQt6 import statements (must be none)."""
    pats = [
        r"^[[:space:]]*from PyQt5",
        r"^[[:space:]]*import PyQt5",
        r"^[[:space:]]*from PySide2",
        r"^[[:space:]]*import PySide2",
        r"^[[:space:]]*from PyQt6",
        r"^[[:space:]]*import PyQt6",
    ]
    args = ["git", "grep", "-nI"]
    for p in pats:
        args += ["-e", p]
    args += ["--", "*.py"]
    r = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    # git grep: rc 0 = matches found (BAD), rc 1 = no match (GOOD)
    hits = [ln for ln in r.stdout.splitlines() if ln.strip()]
    return hits


def check_cubit_toolbar_imports_pyside6() -> list[str]:
    """Each Cubit toolbar module must import PySide6."""
    missing = []
    for rel in CUBIT_TOOLBAR_MODULES:
        p = ROOT / rel
        if not p.is_file():
            missing.append(f"{rel}: missing file")
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        has = "from PySide6" in txt or "import PySide6" in txt
        if not has:
            missing.append(rel)
    return missing


def check_pyproject() -> list[str]:
    """radia pyproject must not depend on PySide6/PyQt or expose old scripts."""
    issues = []
    p = ROOT / "pyproject.toml"
    txt = p.read_text(encoding="utf-8", errors="replace")
    for tok in ('"PySide6', "'PySide6", '"PyQt5', "'PyQt5",
                '"PySide2', "'PySide2", '"PyQt6', "'PyQt6"):
        if tok in txt:
            issues.append(f"pyproject.toml: Qt dependency token appears: {tok}")
    for tok in (
        'radia-ih = "radia.radia_ih:main"',
        'radia-em = "radia.radia_em:main"',
        'radia-pcb = "radia.radia_pcb:main"',
        'radia-streamfunction = "radia.radia_streamfunction:main"',
    ):
        if tok in txt:
            issues.append(f"pyproject.toml: retired script still appears: {tok}")
    return issues


def check_ccm_qt_free() -> tuple[str, list[str]]:
    """cubit_mesh_export.ccm must not import any Qt5 DLL (pefile). Returns (status, issues)."""
    ccm = ROOT / "packages" / "cubit-mesh-export" / "src" / "cubit_mesh_export" / "cubit_mesh_export.ccm"
    if not ccm.is_file():
        return ("skip (no .ccm in package source)", [])
    try:
        import pefile
    except ImportError:
        return ("skip (pefile not installed)", [])
    pe = pefile.PE(str(ccm), fast_load=True)
    pe.parse_data_directories(
        directories=[pefile.DIRECTORY_ENTRY["IMAGE_DIRECTORY_ENTRY_IMPORT"]])
    deps = {e.dll.decode().lower() for e in getattr(pe, "DIRECTORY_ENTRY_IMPORT", [])}
    qt = sorted(d for d in deps if d.startswith("qt"))
    return ("checked", [f".ccm imports Qt DLL(s): {qt}"] if qt else [])


# ----------------------------------------------------------------------
# Headless panel smoke (run in --smoke subprocess on offscreen Qt)
# ----------------------------------------------------------------------
class _StubCubit:
    def get_block_id_list(self):
        return []

    def get_block_name(self, _bid):
        return ""

    def get_volume_count(self):
        return 1


def _run_smoke() -> int:
    """--smoke mode: construct ExportDialog under offscreen Qt."""
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    sys.path.insert(0, str(ROOT / "src" / "radia" / "panels"))
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])

    fails: list[str] = []

    try:
        import radia_export_menu as rem
        stub = _StubCubit()
        for fmt in (rem.FMT_NETGEN, rem.FMT_GMSH, rem.FMT_NASTRAN,
                    rem.FMT_VTK, rem.FMT_FEMEEM, rem.FMT_MEG):
            dlg = rem.ExportDialog(fmt, "C:/temp/model.jou", stub, parent=None)
            cmd = dlg.cubit_command()
            if not (cmd.startswith("export ") and cmd.endswith(" overwrite")):
                fails.append(f"ExportDialog[{fmt}]: bad cmd {cmd!r}")
    except Exception as e:  # noqa: BLE001
        fails.append(f"ExportDialog: {type(e).__name__}: {e}")

    print(json.dumps({"ok": not fails, "fails": fails}))
    return 1 if fails else 0


def check_panel_smoke() -> list[str]:
    """Spawn the toolbar smoke in an isolated offscreen subprocess."""
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    r = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--smoke"],
                       capture_output=True, text=True, env=env, timeout=120)
    line = next((ln for ln in reversed(r.stdout.splitlines()) if ln.strip().startswith("{")), "")
    if not line:
        return [f"toolbar smoke produced no JSON (rc={r.returncode}); stderr tail:\n{r.stderr[-800:]}"]
    data = json.loads(line)
    return list(data.get("fails", [])) if not data.get("ok") else []


# ----------------------------------------------------------------------
def main(argv: list[str]) -> int:
    if "--smoke" in argv:
        return _run_smoke()

    print("=" * 64)
    print("  Cubit toolbar PySide6 / legacy Qt audit")
    print("=" * 64)
    all_issues: list[str] = []

    legacy = check_no_legacy_qt_imports()
    print(f"[{'FAIL' if legacy else 'OK'}] legacy Qt imports (PyQt5/PySide2/PyQt6): "
          f"{len(legacy)} found")
    for h in legacy:
        print("       " + h)
    all_issues += [f"legacy Qt import: {h}" for h in legacy]

    gui_missing = check_cubit_toolbar_imports_pyside6()
    print(f"[{'FAIL' if gui_missing else 'OK'}] Cubit toolbar modules import PySide6")
    all_issues += [f"Cubit toolbar module missing PySide6 import: {m}" for m in gui_missing]

    pyproj = check_pyproject()
    print(f"[{'FAIL' if pyproj else 'OK'}] radia pyproject (no Qt dependency / retired scripts)")
    all_issues += pyproj

    ccm_status, ccm_issues = check_ccm_qt_free()
    print(f"[{'FAIL' if ccm_issues else 'OK'}] cubit_mesh_export.ccm Qt-free  ({ccm_status})")
    all_issues += ccm_issues

    if importlib.util.find_spec("PySide6") is None:
        smoke = []
        print("[OK] headless Cubit toolbar smoke skipped "
              "(PySide6 not installed in this interpreter; Cubit owns the target runtime)")
    else:
        smoke = check_panel_smoke()
        print(f"[{'FAIL' if smoke else 'OK'}] headless Cubit toolbar smoke "
              f"(ExportDialog x6)")
        for s in smoke:
            print("       " + s)
        all_issues += [f"panel smoke: {s}" for s in smoke]

    print("-" * 64)
    if all_issues:
        print(f"RESULT: {len(all_issues)} ISSUE(S) -- Cubit toolbar audit NOT clean")
        return 1
    print("RESULT: CLEAN -- no legacy Qt; Cubit toolbar smoke is healthy/skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
