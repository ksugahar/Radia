---
name: pyside6-health
description: Retired/guard skill for the old PySide6 desktop-panel era. Do NOT uninstall or delete Coreform Cubit's bundled PySide6; Cubit owns that private runtime. For current notebook panels use ipynb-gui-health, and for Cubit deploy use cubit-plugin-install --verify-only plus cubit-smoke-test. Use this only to interpret old PySide6 audit logs or to confirm the Cubit-vs-Radia Python boundary.
---

# pyside6-health

> **Retired / guard-only (2026-06-28).** The Radia operating panels have moved
> to notebook workbenches, and the lab target is now **no PySide6 in the normal
> Radia Python environments** (LAB / 100号機 / mdx / hibino).  The exception is
> **Coreform Cubit's bundled PySide6**, under
> `C:/Program Files/Coreform Cubit*/bin/python3/lib/site-packages`: that belongs
> to Cubit itself and MUST NOT be removed.  Do not run uninstall commands against
> Cubit's embedded Python.  Use `ipynb-gui-health` for notebook panels and
> `cubit-plugin-install --verify-only` + `cubit-smoke-test` for Cubit deploy.

> **Scope narrowed (2026-06-25).** The panel GUIs (radia-ih / em / pcb / motor /
> streamfunction) were **promoted to Jupyter notebook workbenches**
> (`radia.panel_notebook_promotion.v1`); their *paradigm* health now lives in
> the **`ipynb-gui-health`** skill. This skill remains the checker for the parts
> NOT migrated: the **in-Cubit export toolbar** (`radia-export-menu`, manifest
> state `migration-shell` -- `panels/radia_export_menu.py` +
> `register_toolbar.py`, still PySide6 Layer-2) and the **C++ Cubit plugin
> deploy** (`.ccm`/`.pyd`/netgen DLLs, cross-machine). For a notebook panel that
> "does not run / shows no result", use `ipynb-gui-health` instead.

Since radia 4.80.0 the Cubit GUI is **PySide6 (Qt6) only** -- the Qt5 `.ccl`
(RadiaComp.cpp) was deleted; the GUI is the PySide6 toolbar
(`panels/register_toolbar.py` + `panels/radia_export_menu.py`).  Target is
**Coreform Cubit 2025.12**, which bundles PySide6 and cannot load a Qt5 `.ccl`.
See AGENTS.md "Cubit GUI: PySide6-Only -- No Qt5 / PyQt5".

This file is now a historical checklist plus a boundary guard:
1. Normal Radia Python should not require PySide6.
2. Cubit's private bundled PySide6 is protected and is not a cleanup target.
3. Cubit plugin freshness is checked with deploy/smoke commands, not by removing
   Qt packages.

## When to use

- Do not use this as a new release gate.
- Do not use it to justify uninstalling Cubit's embedded PySide6.
- Use it only when reading old PySide6-era notes or confirming that normal
  Python and Cubit's embedded Python are being treated as separate runtimes.

## Layer A -- static + headless (any dev machine with the repo)

```bash
python tools/audit_pyside6_only.py
```

Checks (exit 0 = clean):
1. ZERO real PyQt5/PySide2/PyQt6 import statements in tracked `*.py`.
2. core GUI modules import PySide6.
3. radia pyproject declares PySide6, no PyQt5 dependency.
4. `cubit_mesh_export.ccm` has no Qt5 DLL dependency (pefile import scan).
5. headless panel smoke: IHPanel/EMPanel/PCBPanel construct + ExportDialog
   builds all 6 formats, run in an **isolated offscreen subprocess**.

NOTE: do NOT try to verify panels via `pytest tests/panels/*_qt.py` on LAB --
`tests/conftest.py` does `os.add_dll_directory(mkl_bin)` and a later in-process
PySide6 import crashes with `0xc0000139` (Qt6/MKL DLL clash).  The audit
script's subprocess smoke is the working substitute.  (The pytest suite still
runs on CI where the MKL dir is not injected the same way.)

## Layer B -- LAB Cubit deploy health + end-to-end

The panel plugin (`.ccm` + `.pyd` + netgen DLLs) must match the package source
and the `export` APREPRO command must actually run in Cubit 2025.12.

```bash
# 1. Cubit MUST be closed (it locks .ccm/.pyd). Kill leftovers first:
pwsh -Command "Get-Process -ErrorAction SilentlyContinue | Where-Object { \$_.Name -eq 'coreform_cubit' -or \$_.Name -eq 'cubit' } | Stop-Process -Force; Start-Sleep 2"

# 2. Deploy freshness (sha256 match of every deployed binary vs package):
cubit-plugin-install --verify-only          # expect: every [OK], no mismatch

# 3. If verify-only reports a mismatch -> redeploy (root cause of stale-panel):
cubit-plugin-install --all-users            # cleans old files (incl. stale .ccl), re-verifies

# 4. End-to-end: export actually works in headless Cubit:
cubit-smoke-test                            # expect: [OK] round-trip healthy
```

`cubit-smoke-test` plays a sample `.jou`, runs `export netgen`, and
asserts `source/sink/sibc` + `coil/workpiece/air` labels.  A Cubit `exit=2`
after a successful export is the benign headless-teardown segfault (the .vol is
valid; the test handles it).

Confirm Cubit bundles PySide6 (not PyQt5):

```bash
pwsh -Command "Get-ChildItem 'C:/Program Files/Coreform Cubit*/bin/python3/lib/site-packages' -Directory -Include PySide6,PyQt5 | Select-Object Name"
# expect: PySide6 present, PyQt5 absent
```

## Layer C -- cross-machine (100号機 + mdx via SSH)

100号機 runs the **W: editable + Cubit symlink** install, while mdx runs
the **PyPI** install.  Layer-A static audit runs on source checkouts; the
machine checks below verify deploy state + end-to-end + PySide6 bundle
per host (192.168.11.100 = 100号機; mdx via its ssh alias):

```bash
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
# stop any Cubit (shared lab machine -- coordinate with users first)
Get-Process -ErrorAction SilentlyContinue | Where-Object { $_.Name -eq 'coreform_cubit' } | Stop-Process -Force
Start-Sleep 2
# Cubit ships PySide6, not PyQt5
Get-ChildItem 'C:/Program Files/Coreform Cubit*/bin/python3/lib/site-packages' -Directory -Include PySide6,PyQt5 | Select-Object Name
# plugin deploy sha match + end-to-end
cubit-plugin-install --verify-only
cubit-smoke-test
PS
```

Per AGENTS.md "Distribution Test Policy", a release is not "deploy verified"
until 100号機, mdx, and hibino pass Layer C.  This dovetails with
`python tools/release_qud.py phase9` (cross-machine consistency) and
`done` (definition of done).

## Definition of healthy

- Layer A: `audit_pyside6_only.py` exits 0 (CLEAN).
- Layer B (LAB): `cubit-plugin-install --verify-only` all `[OK]`;
  `cubit-smoke-test` `[OK]`; Cubit has PySide6 / no PyQt5.
- Layer C: same on 100号機 + mdx + hibino.

If any Layer B/C fails on `verify-only` -> the deploy is stale; redeploy with
`cubit-plugin-install --all-users` (Cubit closed) and re-verify.  If Layer A
fails -> a Qt5/PyQt5 reference crept back in; fix it (see AGENTS.md policy) and
do NOT add a `try PySide6 except PyQt5` fallback.

## Related skills / tools

- `tools/audit_pyside6_only.py` -- the Layer-A checker (committed, CI-able).
- `radia-plugin-check` -- broader plugin-freshness triage + 4-layer chain.
- `verify-deploy` -- "are my src/radia edits actually loaded?".
- `panel-cli-diff` / `panel-qt-test` -- panel flag/widget wiring (CI / dev).
- `deploy` / `release-qud` -- the deploy + cross-machine consistency flow.
