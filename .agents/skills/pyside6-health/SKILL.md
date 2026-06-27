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

> **Scope narrowed (2026-06-28).** The panel GUIs (radia-ih / em / pcb /
> motor / streamfunction) were **promoted to Jupyter notebook workbenches**
> (`radia.panel_notebook_promotion.v1`).  Their health now lives in
> **`ipynb-gui-health`**.  This file only preserves the boundary rule for old
> PySide6-era notes: normal Radia Python should not depend on PySide6, while
> Coreform Cubit's private embedded PySide6 must remain intact.

Since radia 4.80.0 the old Qt5 `.ccl` path (RadiaComp.cpp) was removed.  Cubit
2025.12 itself bundles PySide6 for its embedded Python, and the in-Cubit export
toolbar (`panels/register_toolbar.py` + `panels/radia_export_menu.py`) may use
that Cubit-owned runtime.  That does **not** make PySide6 a production
dependency of the normal Radia Python environment.

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

## Layer A -- normal Radia Python boundary

Run this in the normal user Python, not Cubit's embedded Python:

```powershell
python -c "import importlib.util; print('PySide6', 'FOUND' if importlib.util.find_spec('PySide6') else 'MISS'); print('shiboken6', 'FOUND' if importlib.util.find_spec('shiboken6') else 'MISS')"
```

Healthy result: both are `MISS` on LAB / 100号機 / mdx / hibino.  If either is
`FOUND`, remove it from the normal user environment only.  Do not inspect or
modify Coreform Cubit's embedded site-packages in this step.

## Layer B -- Cubit-owned PySide protection + deploy health

The Cubit plugin (`.ccm` + `.pyd` + netgen DLLs) must match the package source
and the `export` APREPRO command must actually run in Cubit 2025.12.  Cubit's
embedded PySide6 is protected and should be present.

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

## Layer C -- cross-machine boundary check

Run the same normal-Python MISS check plus Cubit deploy/protection check per
host (192.168.11.100 = 100号機; mdx / hibino via their ssh aliases):

```bash
cat << 'PS' | ssh 192.168.11.100 'pwsh -ExecutionPolicy Bypass -Command -'
# stop any Cubit (shared lab machine -- coordinate with users first)
Get-Process -ErrorAction SilentlyContinue | Where-Object { $_.Name -eq 'coreform_cubit' } | Stop-Process -Force
Start-Sleep 2
# normal Python should not see PySide6/shiboken6
python -c "import importlib.util; print('PySide6', 'FOUND' if importlib.util.find_spec('PySide6') else 'MISS'); print('shiboken6', 'FOUND' if importlib.util.find_spec('shiboken6') else 'MISS')"
# Cubit ships PySide6, not PyQt5. Do not delete this PySide6.
Get-ChildItem 'C:/Program Files/Coreform Cubit*/bin/python3/lib/site-packages' -Directory -Include PySide6,PyQt5 | Select-Object Name
# plugin deploy sha match + end-to-end
cubit-plugin-install --verify-only
cubit-smoke-test
PS
```

Per AGENTS.md "Distribution Test Policy", a release is not "deploy verified"
until 100号機, mdx, and hibino pass the Cubit deploy checks.  Notebook panel
health is handled by `ipynb-gui-health`.

## Definition of healthy

- Layer A: normal Python reports `PySide6 MISS` and `shiboken6 MISS`.
- Layer B (LAB): `cubit-plugin-install --verify-only` all `[OK]`;
  `cubit-smoke-test` `[OK]`; Cubit has PySide6 / no PyQt5.
- Layer C: same boundary/deploy state on 100号機 + mdx + hibino.

If any Layer B/C fails on `verify-only` -> the deploy is stale; redeploy with
`cubit-plugin-install --all-users` (Cubit closed) and re-verify.  If Layer A
finds PySide6 in normal Python, clean that user Python environment only.  Do
NOT uninstall or delete Cubit's embedded PySide6.

## Related skills / tools

- `ipynb-gui-health` -- current notebook panel health gate.
- `radia-plugin-check` -- broader plugin-freshness triage + deploy chain.
- `deploy` / `release-qud` -- the deploy + cross-machine consistency flow.
