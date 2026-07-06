---
name: verify-deploy
description: Verify that local edits to src/radia/ are actually loaded by notebook workbenches, Cubit toolbar code, and external Python (Python 3.12). Catches stale .pyc, wrong site-packages, UNC vs S: drive confusion, and Cubit cached register_toolbar.py.
---

# verify-deploy

After editing `src/radia/**/*.py`, verify that the changes are
**actually visible** to:

1. The external Python 3.12 (run by notebook workbenches / headless calcs)
2. The Cubit-bundled Python 3.10 (which loads `register_toolbar.py`
   on startup or via Reload Toolbar)
3. The `radia.<app>_notebook` modules and `src/radia/panels/calc_*.py`
   scripts

This skill exists because the LAB editable install has a confusing
twist: `__editable__.radia-4.4.0.pth` points at the **UNC path**
`\\192.168.11.100\work\00_CAE\Radia\01_GitHub\src` (because that is
how the install was performed). The S: drive on the LAB and the W:
drive on 100号機 are both mapped to the same SMB share, so editing
through S: works, but `radia.__file__` reports the UNC form, not S:.

## Phase 1: install state

```bash
python -c "
import radia, sys, os
print('=== install state ===')
print(f'radia.__version__   = {radia.__version__}')
print(f'radia.__file__      = {radia.__file__}')
print(f'radia base dir      = {os.path.dirname(radia.__file__)}')
import importlib.metadata as md
try:
    dist = md.distribution('radia')
    files = dist.files or []
    pth = [f for f in files if str(f).endswith('.pth')]
    print(f'editable .pth       = {pth}')
except Exception as e:
    print(f'metadata lookup     = FAILED: {e}')
"
```

Expected:
- `radia.__file__` is either an S:/W: path or a `\\192.168.11.100\work\` UNC path. Both mean "editable install backed by the SMB-shared repo" — that is OK.
- `radia.__file__` is **NOT** under `Lib/site-packages/radia/` — that would mean a wheel install is shadowing the editable install and edits are silently ignored.

If `radia.__file__` points at site-packages, run:
```bash
pip uninstall -y radia
pip install -e . --no-deps
```

## Phase 2: source-of-truth check (per file)

For every file you just edited, verify that the **mtime + size** of
the module Python is actually loading match what is on disk:

```bash
python -c "
import importlib, os, sys
mods = [
    'radia.bem_coupled_solver',
    'radia.radia_gui_base',
    'radia.radia_ih',
    'radia.bem_inductance',
]
print(f'{\"module\":<35} {\"mtime\":>12} {\"size\":>8}  file')
print('-' * 95)
for name in mods:
    m = importlib.import_module(name)
    p = m.__file__
    try:
        st = os.stat(p)
        print(f'{name:<35} {int(st.st_mtime):>12} {st.st_size:>8}  {p}')
    except OSError as e:
        print(f'{name:<35} ERROR: {e}')
"
```

Compare against the on-disk file with `ls -la src/radia/...`. Times
**must match exactly**. If they do not, the editable install is broken
or `sys.path` has a stale earlier copy of `radia`.

For a calc_*.py panel script (which is loaded by an external Python
subprocess and is **not** a top-level module), check directly:

```bash
ls -la src/radia/panels/calc_inductance.py src/radia/panels/calc_fem_kelvin.py
stat -c "%Y %s %n" src/radia/panels/*.py
```

## Phase 3: stale .pyc check

Python caches compiled .pyc files in `__pycache__/`. Python 3 invalidates
them automatically when the .py mtime is newer, but on SMB shares the
mtime can lag and a stale .pyc may be used briefly.

```bash
# Show the .pyc / .py mtime difference for the files you just edited
python -c "
import os
files = [
    'src/radia/bem_coupled_solver.py',
    'src/radia/radia_gui_base.py',
    'src/radia/radia_ih.py',
    'src/radia/panels/calc_inductance.py',
]
for f in files:
    py_mt = os.path.getmtime(f)
    pyc = os.path.join(os.path.dirname(f), '__pycache__',
                       os.path.basename(f).replace('.py', '.cpython-312.pyc'))
    if os.path.isfile(pyc):
        pyc_mt = os.path.getmtime(pyc)
        delta = pyc_mt - py_mt
        flag = 'OK' if pyc_mt >= py_mt else 'STALE!'
        print(f'  [{flag}] py={int(py_mt)} pyc={int(pyc_mt)} delta={int(delta):+d}s {f}')
    else:
        print(f'  [no .pyc] {f}')
"
```

If anything is `STALE!`, nuke the cache:
```bash
find src/radia -name "__pycache__" -type d -exec rm -rf {} +
```

## Phase 4: signature check (function-level)

The strongest verification: extract the source of a critical function
via `inspect.getsource` and grep for a known string from your edit.
This proves the **loaded function** is the one you wrote, not a stale
binary or shadow copy.

```bash
python -c "
import inspect
from radia.bem_coupled_solver import CoupledBEMSolver, assemble_back_reaction_RHS

# verify the per-DOF f_back fix is in place
src = inspect.getsource(assemble_back_reaction_RHS)
assert 'CoefficientFunction' in src, 'per-DOF f_back missing'
assert 'LinearForm' in src, 'LinearForm assembly missing'
assert 'wp_J[j, k]' in src, 'wp_J indexing missing'
print('OK: assemble_back_reaction_RHS contains per-DOF assembly')

# verify the iteration loop uses both J_re AND J_im (sign-correct)
src = inspect.getsource(CoupledBEMSolver.solve)
assert 'J_re' in src and 'J_im' in src, 'real+imag split missing'
assert 'L_air' in src, 'L_air baseline missing'
print('OK: CoupledBEMSolver.solve uses real+imag back-reaction')
"
```

Adapt the asserts for whatever you just changed. The point is: if the
assertion fires, you know the import is reading a different file than
the one you edited.

## Phase 5: end-to-end CLI smoke test

Run the actual entry point you care about and check that the new
output strings appear:

```bash
# IH panel BEM coupled, headline numbers
python src/radia/panels/calc_inductance.py \
    --vol src/radia/panels/samples/radia_model.vol \
    --frequency 50000 --coil-sigma 5.8e7 \
    --workpiece workpiece --impedance-model dowell \
    --sigma 5.8e7 --half-thickness 0.005 --mu-r 1 \
    --msh-output /tmp/_verify.msh 2>&1 | tail -1 | python -c "
import sys, json
d = json.loads(sys.stdin.read())
assert 'coupled_dL_H' in d, 'coupled solver did not run'
dL = d['coupled_dL_H'] * 1e9
print(f'BEM coupled dL = {dL:+.3f} nH (expect copper Lenz < 0)')
assert dL < 0, f'wrong sign for copper: {dL:+.3f} nH'
print('OK')
"
```

If any of these assertions fail, the editable install is not picking
up your edits. Go back to Phase 1.

## Phase 6: GUI display check (radia_gui_base.py)

`radia_gui_base.py` formats the result dict into the panel output
window. The formatter changes are only seen when an external Python
process imports `radia.radia_gui_base`. Easy check:

```bash
python -c "
import inspect
from radia.radia_gui_base import AnalysisWindow
src = inspect.getsource(AnalysisWindow._on_finished)
# verify the new BEM coupled display strings are in place
assert 'L (air)' in src, 'BEM coupled display block missing'
assert 'delta L' in src, 'delta L display missing'
print('OK: AnalysisWindow._on_finished has the new display block')
"
```

## Phase 7: Cubit-side reload check

The Cubit-bundled Python 3.10 loads `register_toolbar.py` only at
startup (or via `Solve > Reload Panels`). After editing
`register_toolbar.py` you MUST either:

1. Restart Cubit completely (`taskkill //F //IM coreform_cubit.exe`
   then relaunch), or
2. Click `Solve > Reload Panels` from the menu

The panel debug log records the load timestamp:

```bash
head -5 C:/radia_panel_log.txt
# The timestamp on the first line should match your Cubit launch /
# reload time. If it does not, Cubit is still running the old module.
```

## Phase 7b: Unified panel debug log check

All Radia GUI components write to **one** file:
`C:/radia_panel_log.txt` (or `$TMPDIR/radia_panel_log.txt` on
non-Windows). Each line is tagged with **timestamp + user@host +
source component** so you can read a single Cubit session as one
continuous log AND tell which user on which machine produced it:

```
[2026-04-12 14:30:12.345] [ksugahar@LAB         ] (cubit       ) register_toolbar.py loaded
[2026-04-12 14:30:18.892] [ksugahar@LAB         ] (cubit       ) _launch_radia_ngsolve: ENTER
[2026-04-12 14:30:25.103] [ksugahar@LAB         ] (ih-window   ) _on_run: cmd=...
[2026-04-12 14:30:25.567] [ksugahar@LAB         ] (inductance  ) calc_main: args={...}
[2026-04-12 14:30:37.842] [ksugahar@LAB         ] (inductance  ) calc_main: done 12.3s keys=[...]
[2026-04-12 14:30:37.901] [ksugahar@LAB         ] (ih-window   )   L = 87.81 nH
```

When the student (Kubota) reports a problem from 100号機, his lines
will show ``[kubota@KUBOTA-PC      ]`` so the LAB and 100号機 sessions
are immediately distinguishable in a single shared log dump.

The user@host tag is captured **once at process start** (in
``init_panel_log``) so it does not change mid-session even if the
environment shifts. The session start banner also records
``user=...`` and ``host=...`` in plain text for grep-ability.

Components and where they run:

| Tag           | Process                | Source file                  |
|---------------|------------------------|------------------------------|
| `cubit`       | Cubit Python 3.10      | `panels/register_toolbar.py` |
| `ih-window`   | external Python 3.12   | `radia_gui_base.py`          |
| `inductance`  | external Python 3.12   | `panels/calc_inductance.py`  |
| `fem_kelvin`  | external Python 3.12   | `panels/calc_fem_kelvin.py`  |
| `heating_bem` | external Python 3.12   | `panels/calc_heating_bem.py` |
| (others)      | external Python 3.12   | other `panels/calc_*.py`     |

The shared writer is `src/radia/panels/panel_log.py`. Truncation is
done **only** by `register_toolbar.py` on Cubit session start; all
other processes append. So one Cubit session = one continuous log
file across all subprocesses.

Verify:

```bash
# Quick smoke check that all 3 components can write
python -c "
import sys, os
sys.path.insert(0, 'src/radia/panels')
from panel_log import init_panel_log, panel_log, PANEL_LOG_PATH
print(f'log path: {PANEL_LOG_PATH}')
init_panel_log('verify', truncate=False, banner=False)
panel_log('verify-deploy phase 7b smoke check')
print('OK: panel_log writable')
"

# Tail the last session to see what happened
tail -50 C:/radia_panel_log.txt 2>/dev/null || tail -50 /tmp/radia_panel_log.txt
```

If `C:/radia_panel_log.txt` cannot be written (permission denied),
the panels still run — log writes are best-effort and never raise.
But you lose the post-mortem record. Check that the user running
Cubit has write permission on `C:/`.

## Phase 8: 100号機 side reflection (NAS-shared repo)

Because `S:\Radia\01_GitHub` (LAB) and `W:\00_CAE\Radia\01_GitHub`
(100号機 NAS mount) are the **same SMB share**, edits made on LAB are
**immediately visible** on 100号機. There is NO need for a separate
SMB copy step. However:

- The 100号機 Cubit holds its own register_toolbar.py cache. After
  editing, do `ssh 100 'Stop-Process -Name coreform_cubit -Force'`.
- If you also installed a wheel into 100号機's site-packages
  (`pip install --force-reinstall ...`), THAT will shadow the
  editable install. Uninstall via `ssh 100 'pip uninstall -y radia'`
  and rely on the editable install via the SMB share.

## Quick one-liner summary

```bash
# After every nontrivial edit, run this BEFORE asking the user to test:
python -c "
import os, importlib
for name in ['radia.bem_coupled_solver','radia.radia_gui_base','radia.radia_ih']:
    m = importlib.import_module(name)
    print(f'  {name:<32} mtime={int(os.path.getmtime(m.__file__))} {m.__file__}')
" && find src/radia -name "__pycache__" -type d | head -5
```

## Arguments

| Argument | Action |
|----------|--------|
| (empty) | Run phases 1-6 (no Cubit interaction) |
| `cubit` | Phases 1-6 + Phase 7 (asks user to restart Cubit) |
| `100` | Phases 1-6 + Phase 8 (verify 100号機) |
| `clean` | Nuke all __pycache__ + run phases 1-6 |
