#!/usr/bin/env python3
"""Audit: Cubit GUI has no legacy Qt and optional PySide6 panels stay healthy.

Enforces the Radia policy "no Qt5 / PyQt5 fallbacks" while allowing release
targets that intentionally do not install PySide6.  When PySide6 is installed,
the legacy desktop panels must still construct headlessly.  When PySide6 is not
installed, the smoke is skipped and the release remains valid because notebook
panels and the Cubit export plugin do not require PySide6.

Static checks:
  1. ZERO real PyQt5 / PySide2 / PyQt6 import statements in tracked *.py
     (lint-rule regexes and "no PyQt5 fallback" comments do NOT count --
     only actual ``import``/``from ... import`` lines).
  2. Core GUI modules import PySide6.
  3. radia pyproject.toml declares PySide6 and lists no PyQt5 dependency.
  4. cubit_mesh_export.ccm (if present in the package source) has NO Qt5 DLL
     dependency (pefile import-table scan).

Headless check (isolated subprocess, offscreen Qt):
  5. If PySide6 is installed, IHPanel / EMPanel / PCBPanel construct under
     PySide6 offscreen, and ExportDialog builds all 6 formats emitting valid
     ``export ... overwrite`` commands.  If PySide6 is absent, this is a
     release-policy skip, not a failure.

Run:  python tools/audit_pyside6_only.py
Exit 0 = no legacy Qt and optional panel surface healthy/skipped; non-zero =
issues (listed).

The machine/Cubit half (cubit-plugin-install --verify-only, cubit-smoke-test,
cross-machine SSH to 100号機 / mdx) lives in the `pyside6-health` skill, which
calls this script for the static/headless portion.
"""
from __future__ import annotations

import json
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

GUI_MODULES = [
    "src/radia/radia_gui_base.py",
    "src/radia/radia_ih.py",
    "src/radia/radia_em.py",
    "src/radia/radia_motor.py",
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


def check_gui_imports_pyside6() -> list[str]:
    """Each core GUI module must import PySide6 (directly or be a base user)."""
    missing = []
    # radia_pcb inherits Qt from radia_gui_base, so only the base + the
    # in-Cubit toolbars + windows need a direct PySide6 import.
    must_import = {
        "src/radia/radia_gui_base.py",
        "src/radia/panels/radia_export_menu.py",
        "src/radia/panels/register_toolbar.py",
    }
    for rel in GUI_MODULES:
        p = ROOT / rel
        if not p.is_file():
            continue
        txt = p.read_text(encoding="utf-8", errors="replace")
        has = "from PySide6" in txt or "import PySide6" in txt
        if rel in must_import and not has:
            missing.append(rel)
    return missing


def check_pyproject() -> list[str]:
    """radia pyproject must declare PySide6 and not depend on PyQt5."""
    issues = []
    p = ROOT / "pyproject.toml"
    txt = p.read_text(encoding="utf-8", errors="replace")
    # Crude but robust: PySide6 must appear as a quoted dependency token.
    if "PySide6" not in txt:
        issues.append("pyproject.toml: PySide6 not declared")
    # PyQt5 as a dependency token (quoted) is forbidden.  A bare mention in
    # a comment is allowed, so look for a quoted requirement form.
    for tok in ('"PyQt5', "'PyQt5"):
        if tok in txt:
            issues.append("pyproject.toml: PyQt5 appears as a dependency")
            break
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
    """--smoke mode: construct panels + ExportDialog under offscreen Qt."""
    os.environ["QT_QPA_PLATFORM"] = "offscreen"
    sys.path.insert(0, str(ROOT / "src" / "radia"))
    sys.path.insert(0, str(ROOT / "src" / "radia" / "panels"))
    from PySide6.QtWidgets import QApplication
    QApplication.instance() or QApplication([])

    fails: list[str] = []

    for modname, clsname in (("radia_ih", "IHPanel"),
                             ("radia_em", "EMPanel"),
                             ("radia_pcb", "PCBPanel")):
        try:
            mod = __import__(modname)
            getattr(mod, clsname)()  # construct
        except Exception as e:  # noqa: BLE001
            fails.append(f"{clsname} construct: {type(e).__name__}: {e}")

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
    """Spawn the smoke in an isolated offscreen subprocess (no conftest DLL pollution)."""
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    r = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--smoke"],
                       capture_output=True, text=True, env=env, timeout=120)
    line = next((ln for ln in reversed(r.stdout.splitlines()) if ln.strip().startswith("{")), "")
    if not line:
        return [f"panel smoke produced no JSON (rc={r.returncode}); stderr tail:\n{r.stderr[-800:]}"]
    data = json.loads(line)
    return list(data.get("fails", [])) if not data.get("ok") else []


# ----------------------------------------------------------------------
def main(argv: list[str]) -> int:
    if "--smoke" in argv:
        return _run_smoke()

    print("=" * 64)
    print("  Optional PySide6 / legacy Qt audit")
    print("=" * 64)
    all_issues: list[str] = []

    legacy = check_no_legacy_qt_imports()
    print(f"[{'FAIL' if legacy else 'OK'}] legacy Qt imports (PyQt5/PySide2/PyQt6): "
          f"{len(legacy)} found")
    for h in legacy:
        print("       " + h)
    all_issues += [f"legacy Qt import: {h}" for h in legacy]

    gui_missing = check_gui_imports_pyside6()
    print(f"[{'FAIL' if gui_missing else 'OK'}] core GUI modules import PySide6")
    all_issues += [f"GUI module missing PySide6 import: {m}" for m in gui_missing]

    pyproj = check_pyproject()
    print(f"[{'FAIL' if pyproj else 'OK'}] radia pyproject (PySide6 declared, no PyQt5 dep)")
    all_issues += pyproj

    ccm_status, ccm_issues = check_ccm_qt_free()
    print(f"[{'FAIL' if ccm_issues else 'OK'}] cubit_mesh_export.ccm Qt-free  ({ccm_status})")
    all_issues += ccm_issues

    if importlib.util.find_spec("PySide6") is None:
        smoke = []
        print("[OK] headless panel smoke skipped "
              "(PySide6 not installed; release policy permits no-PySide6 targets)")
    else:
        smoke = check_panel_smoke()
        print(f"[{'FAIL' if smoke else 'OK'}] headless panel smoke "
              f"(IHPanel/EMPanel/PCBPanel construct; ExportDialog x6)")
        for s in smoke:
            print("       " + s)
        all_issues += [f"panel smoke: {s}" for s in smoke]

    print("-" * 64)
    if all_issues:
        print(f"RESULT: {len(all_issues)} ISSUE(S) -- optional PySide6 audit NOT clean")
        return 1
    print("RESULT: CLEAN -- no legacy Qt; optional PySide6 panel smoke is healthy/skipped")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
