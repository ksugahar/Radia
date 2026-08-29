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

Interface checks:
  5. The official WorkflowToolbar remains complete and startup registers the
     complementary Export menu through emclaro, never through a PySide6 QMenu.

Headless check (isolated subprocess, offscreen Qt):
  6. If PySide6 is installed, ExportDialog builds all 6 formats emitting valid
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
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
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


def check_deployed_panel_source(
        cubit_file: Path | None = None,
        expected_register: Path | None = None) -> tuple[str, list[str]]:
    """Reject a registered startup shim that loads another checkout.

    The package-level verifier proves that *its own* deployment is internally
    consistent.  A developer audit also has to prove that Cubit loads the
    checkout being audited; otherwise an old release worktree can pass while
    the current menu fix is never executed.
    """
    cubit_file = cubit_file or (Path.home() / ".cubit")
    expected_register = expected_register or (
        ROOT / "src" / "radia" / "panels" / "register_toolbar.py")
    if not cubit_file.is_file():
        return ("skip (no ~/.cubit)", [])

    cubit_text = cubit_file.read_text(encoding="utf-8", errors="replace")
    begin = cubit_text.find("## BEGIN radia toolbar")
    end = cubit_text.find("## END radia toolbar")
    if begin < 0 or end < begin:
        return ("skip (Radia startup not registered)", [])
    match = re.search(r'play\s+"([^"]+)"', cubit_text[begin:end])
    if not match:
        return ("checked", [f"invalid Radia startup block: {cubit_file}"])

    startup = Path(match.group(1).replace("/", os.sep))
    if not startup.is_file():
        return ("checked", [f"registered startup shim is missing: {startup}"])
    startup_text = startup.read_text(encoding="utf-8", errors="replace")
    expected = str(expected_register.resolve()).replace("\\", "/").lower()
    actual_text = startup_text.replace("\\", "/").lower()
    if expected not in actual_text:
        loaded = re.findall(r"[^'\"]+register_toolbar\.py", actual_text)
        loaded_text = loaded[-1] if loaded else "unknown"
        return ("checked", [
            "Cubit startup loads a different checkout: "
            f"expected={expected}; loaded={loaded_text}"
        ])
    return ("checked", [])


def check_official_toolbar_contract() -> list[str]:
    """Validate the WorkflowToolbar and Claro-owned menu contracts."""
    toolbar_root = ROOT / "src" / "radia" / "panels" / "cubit_toolbar"
    template = toolbar_root / "toolbars" / "radia_export_toolbar.ttb.tmpl"
    register = ROOT / "src" / "radia" / "panels" / "register_toolbar.py"
    export_menu = ROOT / "src" / "radia" / "panels" / "radia_export_menu.py"
    issues: list[str] = []
    if not template.is_file():
        return [f"missing WorkflowToolbar template: {template}"]
    try:
        root = ET.fromstring(template.read_text(encoding="utf-8"))
    except (OSError, ET.ParseError) as exc:
        return [f"invalid WorkflowToolbar template: {exc}"]
    if root.tag != "WorkflowToolbar" or root.attrib.get("name") != "Radia Export":
        issues.append("toolbar root must be WorkflowToolbar name='Radia Export'")
    buttons = root.findall("WTButton")
    if len(buttons) != 6:
        issues.append(f"expected 6 official toolbar buttons, found {len(buttons)}")

    for filename in root.findall(".//filename"):
        relative = (filename.text or "").split("@TOOLBAR_INSTALL_DIR@/", 1)[-1]
        if not relative or not (toolbar_root / relative).is_file():
            issues.append(f"missing toolbar launcher: {relative!r}")

    register_text = register.read_text(encoding="utf-8", errors="replace")
    menu_text = export_menu.read_text(encoding="utf-8", errors="replace")
    if "radia_export_menu.install_menu()" not in register_text:
        issues.append("~/.cubit startup does not install the Claro Export menu")
    if "emclaro.add_to_menu" not in menu_text:
        issues.append("Export menu does not use Cubit's official Claro API")
    if "QMenu(" in menu_text or "menu_bar.addMenu" in menu_text:
        issues.append("Export menu still injects an unsupported PySide6 QMenu")
    return issues


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


def _find_cubit_python() -> Path | None:
    """Return Cubit's embedded Python, which owns production PySide6."""
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
        program_files = Path(os.environ.get("ProgramFiles", r"C:\Program Files"))
        candidates = sorted(
            program_files.glob("Coreform Cubit */bin/python3/python.exe"),
            reverse=True,
        )
        if candidates:
            return candidates[0]
    return None


def _qt_runtime_python() -> tuple[Path | None, str]:
    if importlib.util.find_spec("PySide6") is not None:
        return Path(sys.executable), "current Python"
    cubit_python = _find_cubit_python()
    if cubit_python is not None:
        return cubit_python, "Cubit embedded Python"
    return None, "unavailable"


def check_panel_smoke(python_executable: Path) -> list[str]:
    """Spawn the toolbar smoke in the selected PySide6 runtime."""
    env = os.environ.copy()
    env["QT_QPA_PLATFORM"] = "offscreen"
    r = subprocess.run([str(python_executable), str(Path(__file__).resolve()), "--smoke"],
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

    deploy_status, deploy_issues = check_deployed_panel_source()
    print(f"[{'FAIL' if deploy_issues else 'OK'}] Cubit startup source "
          f"matches audited checkout  ({deploy_status})")
    for issue in deploy_issues:
        print("       " + issue)
    all_issues += deploy_issues

    toolbar_contract = check_official_toolbar_contract()
    print(f"[{'FAIL' if toolbar_contract else 'OK'}] official Coreform "
          "WorkflowToolbar + Claro Export menu contract")
    for issue in toolbar_contract:
        print("       " + issue)
    all_issues += toolbar_contract

    qt_python, qt_runtime = _qt_runtime_python()
    if qt_python is None:
        smoke = []
        print("[OK] headless Cubit toolbar smoke skipped "
              "(no PySide6 runtime or Cubit embedded Python found)")
    else:
        smoke = check_panel_smoke(qt_python)
        print(f"[{'FAIL' if smoke else 'OK'}] headless Cubit toolbar smoke "
              f"(ExportDialog x6; {qt_runtime}: {qt_python})")
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
