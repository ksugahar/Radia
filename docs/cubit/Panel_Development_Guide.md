# Cubit Panel Development Guide

How to add a new solver panel to the Radia Cubit Solve menu.

## Menu Structure

```
Solve (Python, register_toolbar.py):
  Radia-NGSolve...     -> export .vol, launch standalone app
  Generate Coil...     -> CoilBuilder script -> STEP -> import
  --------
  Reload Panels        -> re-read register_toolbar.py (debug)
```

## Important Rules

- **ASCII only** in all .py files loaded by Cubit (cp932 Japanese Windows)
- **Import QMenu** if using submenus (PySide6.QtWidgets; Cubit 2025.12 ships PySide6)
- **No Qt in calc_*.py** -- subprocess scripts must not import PySide6 or any Qt binding
- **No cubit import in calc_*.py** -- .vol file is the sole interface
- `cubit-plugin-install` generates the Cubit startup shim under
  `%ProgramData%/Radia/Cubit/` for `--all-users` installs (or
  `%LOCALAPPDATA%/Radia/Cubit/` for current-user installs). Do not edit
  `src/radia/panels/startup.py` by hand.
- `cubit-plugin-install --verify-only --all-users` must pass before a
  first install is considered complete; it checks both plugin hashes and
  panel startup registration.

## Architecture

```
register_toolbar.py (Cubit GUI Python, Qt)
  |
  |-- YourDialog(QDialog)
  |     1. Display UI (Qt widgets)
  |     2. cubit.cmd('export netgen "model.vol" order 2 overwrite')
  |     3. QProcess -> [external_python, calc_your.py, --vol model.vol, ...]
  |     4. Parse JSON result
  |     5. Display result in dialog
  |
  v
calc_your.py (External Python, NO Qt, NO cubit)
  |-- from calc_common import setup_paths, calc_main, progress, MU_0
  |-- def solve_your(vol_file, ...): ... return dict
  |-- calc_main(solve_your, parser)
```

**Key rules**:
- `calc_*.py` runs in external Python (with NGSolve). It must NOT import Qt.
- `calc_*.py` must NOT import cubit. The `.vol` file is the sole interface.
- Mesh curving is done at export time (`export netgen ... order N`).

## Step-by-Step: Adding a New Panel

### 1. Create `calc_your.py`

```python
"""
Your solver description.

Usage:
    python calc_your.py --vol model.vol --param1 value

Outputs JSON to stdout.
"""

import argparse
import os
import sys

# Shared utilities
_this_dir = os.path.dirname(os.path.abspath(__file__))
if _this_dir not in sys.path:
    sys.path.insert(0, _this_dir)

from calc_common import (setup_paths, progress, calc_main, MU_0, NU_0)


def solve_your(vol_file, param1=1.0):
    """Your solver.

    Args:
        vol_file: Netgen .vol file with material/boundary labels.

    Returns:
        dict with result keys
    """
    from ngsolve import Mesh, Integrate, CF, BND

    setup_paths()

    # 1. Load mesh from .vol (labels embedded by export netgen)
    mesh = Mesh(vol_file)
    progress("MESH", f"{mesh.GetNE()} elements")

    # 2. Solve
    progress("SOLVE", "starting...")
    result_value = 42.0  # your computation

    # 3. Return JSON-serializable dict
    return {
        "result": result_value,
        "n_elements": mesh.GetNE(),
    }


def main():
    parser = argparse.ArgumentParser(description="Your solver")
    parser.add_argument("--vol", required=True, help="Netgen .vol file")
    parser.add_argument("--param1", type=float, default=1.0, help="Parameter")
    parser.add_argument("--output", default="", help="JSON output file")

    def run(args):
        return solve_your(args.vol, args.param1)

    calc_main(run, parser)


if __name__ == "__main__":
    main()
```

### 2. Create Dialog in `register_toolbar.py`

```python
class YourDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Your Panel")
        self._ext_python = _find_external_python()
        self._process = None
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        # ... add widgets ...

        self.solve_btn = QPushButton("Solve")
        self.solve_btn.clicked.connect(self._solve)
        layout.addWidget(self.solve_btn)

        self.result_label = QLabel("")
        layout.addWidget(self.result_label)

    def _solve(self):
        if not self._ext_python:
            QMessageBox.warning(self, "Error", "External Python not found.")
            return

        self.solve_btn.setEnabled(False)
        self.solve_btn.setText("Solving...")

        # Export mesh to .vol (curving done here, at export time)
        tmpdir = tempfile.mkdtemp(prefix="radia_your_")
        vol_file = os.path.join(tmpdir, "model.vol").replace("\\", "/")
        cubit.cmd(f'export netgen "{vol_file}" order 2 overwrite')

        # JSON output file
        self._result_json = os.path.join(tmpdir, "result.json").replace("\\", "/")

        calc_script = os.path.join(_this_dir, "calc_your.py")
        args = [
            calc_script,
            "--vol", vol_file,
            "--output", self._result_json,
        ]

        self._process = QProcess(self)
        self._process.finished.connect(self._on_finished)
        self._process.start(self._ext_python, args)

    def _on_finished(self, exit_code, exit_status):
        self.solve_btn.setEnabled(True)
        self.solve_btn.setText("Solve")
        self._process = None

        data = None
        if os.path.exists(self._result_json):
            try:
                with open(self._result_json, "r") as f:
                    data = json.load(f)
            except Exception:
                pass

        if data is None:
            self.result_label.setText(f"Error (exit code {exit_code})")
            return
        if "error" in data:
            self.result_label.setText(f"Error: {data['error']}")
            return

        self.result_label.setText(f"Result: {data['result']}")
```

### 3. Register in Menu

In `register_menu()`:

```python
action_your = QAction("Your Panel...", main_window)
action_your.setStatusTip("Description")
def _show_your():
    if not hasattr(main_window, '_radia_your_dlg') or main_window._radia_your_dlg is None:
        main_window._radia_your_dlg = YourDialog(main_window)
    main_window._radia_your_dlg.show()
    main_window._radia_your_dlg.raise_()
action_your.triggered.connect(_show_your)
radia_menu.addAction(action_your)
```

Add cleanup for stale dialogs in the reload handler:

```python
for attr in ('_radia_ind_dlg', '_radia_fem_dlg', '_radia_your_dlg'):
    ...
```

## Subprocess Output Protocol

### stdout: JSON result (last line)

`calc_main()` handles this automatically. The solver function returns a dict,
which is serialized as JSON on stdout and optionally to `--output` file.

### stderr: Progress messages

Use `progress(TAG, msg)` from `calc_common`:

```python
progress("MESH", "42 elements exported")    # -> stderr: MESH:42 elements exported
progress("SOLVE", "iteration 3, tol=1e-5")  # -> stderr: SOLVE:iteration 3, tol=1e-5
```

Tags are parsed in `_on_stderr()` handlers in register_toolbar.py.

### Error handling

`calc_main()` catches exceptions and returns `{"error": str(e)}`.
Full traceback is written to stderr for debugging.

## Block Naming Conventions

Cubit blocks map to NGSolve boundary/material labels. Use these standard names:

### Volume Blocks (materials)

| Block Name | Purpose | Required |
|-----------|---------|----------|
| `coil` | Coil conductor | Yes (IH) |
| `workpiece` | Workpiece (SIBC/ESIM target) | Optional |
| `air` | Air domain | Yes (FEM) |
| `kelvin` | Kelvin transform domain (exterior) | Optional (Add Kelvin button) |
| `yoke` | Iron yoke (accelerator magnets) | Optional |

### Surface Blocks (boundaries)

| Block Name | Purpose | Required |
|-----------|---------|----------|
| `source` | Current injection face (T0=1) | Optional* |
| `sink` | Current return face (T0=0) | Optional* |
| `wp_surface` | Workpiece-air interface (Robin BC) | Auto-created |
| `kelvin_int` | Interior Kelvin sphere surface | Auto-created |
| `kelvin_ext` | Exterior Kelvin sphere surface | Auto-created |
| `outer` | Outer boundary (far field) | Auto-created |
| `GND` | Dirichlet at infinity (vertex at Kelvin center) | Optional (HCurl), Essential (H1) |

*source/sink: Required for T0 technique. Without them, J_theta fallback is used (axisymmetric torus only).

### Auto-Created Blocks

The panel dialog auto-creates these blocks when "Run Journal" or "Add Kelvin" is clicked:

- **wp_surface**: Shared faces between workpiece and air volumes
- **kelvin_int/kelvin_ext**: Interior/exterior Kelvin sphere surfaces
- **outer**: Free surfaces on air or kelvin volumes
- **GND**: Vertex at exterior Kelvin sphere center (maps to physical infinity).
  Essential for H1 (scalar potential), optional for HCurl (gauge reg suffices).

### Naming Rules

1. Use **lowercase** names (Cubit is case-insensitive, but NGSolve preserves case)
2. Use **underscores** for multi-word names (`wp_surface`, not `WP-Surface`)
3. Volume blocks contain volumes; surface blocks contain surfaces/tris
4. One block per physical region (no overlapping blocks for the same purpose)
5. Block names become NGSolve `mesh.GetMaterials()` (volumes) or `mesh.GetBoundaries()` (surfaces)

### Example Journal (IH with Kelvin)

```
# Geometry
create cylinder radius 0.003 height 0.2
create cylinder radius 0.01 height 0.04

# Mesh
volume 1 scheme tetmesh
volume 1 size 0.003
mesh volume 1

volume 2 scheme tetmesh
volume 2 size 0.002
mesh volume 2

# Blocks
block 1 add volume 1
block 1 name "coil"

block 2 add volume 2
block 2 name "workpiece"

# Air sphere (will be webcut by Add Kelvin)
create sphere radius 0.06
subtract volume 1 2 from volume 3 keep
volume 3 scheme tetmesh
volume 3 size 0.01
mesh volume 3
block 3 add volume 3
block 3 name "air"

# Source/sink on coil ends
# (faces at z=0 and z=0.2)
block 4 add surface 1
block 4 name "source"
block 5 add surface 2
block 5 name "sink"
```

Then click "Add Kelvin" in the panel to auto-add the Kelvin sphere pair.

## Files

| File | Runs in | Purpose |
|------|---------|---------|
| `src/radia/panels/register_toolbar.py` | Cubit GUI Python | Qt dialogs + menu |
| `src/radia/panels/calc_common.py` | External Python | Shared utilities |
| `src/radia/panels/calc_volume.py` | External Python | Volume integration |
| `src/radia/panels/calc_surface.py` | External Python | Surface area |
| `src/radia/panels/calc_inductance.py` | External Python | BEM inductance |
| `src/radia/panels/calc_fem_kelvin.py` | External Python | FEM-SIBC + Kelvin |
| generated `radia_startup.py` | Cubit GUI Python | Startup shim generated by `cubit-plugin-install` |
| `src/radia/panels/startup.py` | Cubit GUI Python | Legacy compatibility shim only |
| `install_panels.py` | System Python | Installer |
