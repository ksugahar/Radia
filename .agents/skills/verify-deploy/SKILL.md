---
name: verify-deploy
description: Verify that Radia edits are loaded from the intended editable checkout by external Python and Simulink, while Cubit loads only its deployed embedded-Python toolbar files.
---

# Verify Deploy

Use this after changing `src/radia`, Simulink application code, or the Cubit
toolbar. The supported runtime boundary is:

- external Python 3.12: editable `radia` and headless application code;
- MATLAB/Simulink: the same external source through the checked MATLAB/MEX
  entry point and tracked Radia library;
- Coreform Cubit: deployed `.ccm`, curver module, startup hook, and toolbar in
  Cubit's private Python/PySide6 runtime.

The retired desktop PySide applications and notebook workbenches are not a
deployment target.

## 1. Identify the loaded Python source

Run from the machine being verified:

```powershell
python -c "import pathlib,radia; print(radia.__version__); print(pathlib.Path(radia.__file__).resolve())"
python -c "import pathlib,radia.simulink.application as m; print(pathlib.Path(m.__file__).resolve())"
```

Both paths must resolve to the source recorded as the editable intent for this
host and interpreter (`python tools/release_quad.py verify-editable`, or
`python tools/verify_lab_editable.py` on LAB), not a wheel under
`site-packages`. LAB and 100号機 use editable installs. A pointer without a
record is UNVERIFIED, not wrong: record it if it is intended
(`repoint --record-current --reason "<why>"`). If a different source is
intended, move it explicitly and let the record follow:

```powershell
python tools/release_quad.py repoint --package radia --source <checkout> --reason "<why>"
```

`repoint` records the previous pointer, runs `pip install -e`, verifies a
fresh-process import and records the new intent. It does not uninstall first
and does not stop processes; a locked entry point is reported for you to
resolve at a quiet boundary. Reconnect long-lived MCP and MATLAB processes
after changing the editable source. A process keeps imported modules and
registered tool objects until it reloads or restarts.

## 2. Verify the Simulink application boundary

Run the repository's structural and lifecycle checks:

```powershell
python tools/audit_new_panel_contract.py
python -m pytest tests/test_simulink_application.py -q
```

Then use `simulink-app-health`. Open the exact tracked `.slx`, confirm its
resolved `FileName`, run `model_read`, `model_check`, and the focused MATLAB
test. Generated scratch models do not prove that the tracked production model
loads the edited source.

For a solver-backed run, inspect `run.log` and `result.json`. Provenance must
name the expected checkout/backend, and spatial runs must index their checked
`.vol` report and GMSH post-processing artifact.

## 3. Verify the Cubit boundary

Do not import Cubit's PySide6 into external Radia Python. Verify and install the
toolbar with the package-owned commands:

```powershell
check-vol --help
cubit-plugin-install
python tools/audit_pyside6_only.py
```

Inside Cubit, use **Reload Toolbar**, then exercise one export action. Confirm
that the deployed startup hook, toolbar module, `.ccm`, and curver binary have
the expected version/hash. Use `radia-plugin-check` for the complete deployed
inventory and a real export/check cycle.

## 4. Cross-machine acceptance

- LAB and 100号機: canonical editable source and fast application checks.
- mdx: isolated CI/compute environment; never infer its state from LAB files.
- hibino: optional long optimization/validation compute host.
- Release verification: use `release-quad`; every result must be attributable to
  one immutable commit SHA.

Reject the deployment when any machine resolves a different source tree, a
long-lived process exposes a stale tool list, the tracked Simulink model cannot
be reopened cleanly, or Cubit binaries and toolbar files disagree.
