"""
Radia EM (Electromagnet) analysis window -- generator-driven (POLICY
2026-05-30 in radia_gui_base.py).

Architecture:
    EMPanel = wp_vol + Method combo + QStackedWidget over 5 sub-panels.
    Each sub-panel is a ModePanel that calls bind_argparser() with the
    matching calc_*.py.  Method switch reveals the corresponding
    sub-panel; per-method widget state survives across switches via
    save_state's method-prefixed namespace.

    Omega / A-Phi   -> _AccelMagnetPanel(calc_accel_magnet.py)
    MSC             -> _MSCPanel       (calc_accel_msc.py)
    Kelvin Benchmark-> _KelvinBenchPanel(calc_kelvin_benchmark.py)
    Clebsch hodograph -> _ClebschPanel(calc_clebsch_hodograph.py)

The old union-of-widgets EMPanel with method-driven visibility toggles
(437 lines as of 2026-05-30) was replaced by this composite design.
Each sub-panel only knows its own argparse + its own dispatch, so:
  - adding `--new-arg` to one calc surfaces a widget only in THAT
    sub-panel (no other panel touched)
  - cross-method widget leakage is impossible by construction
  - panel-cli-diff drift surface is zero per sub-panel

Usage:
    python -m radia.radia_em model.vol
    python radia_em.py model.vol
"""

import sys
import os

TITLE = "Electromagnet"
REQUIRED_LABELS = []
OPTIONAL_LABELS = ["iron", "magnetic", "air", "kelvin", "sphere",
                   "kelvin_int", "kelvin_ext"]
OPTIONAL_FILES = {"Coil script": "Python (*.py)",
                  "Hysteresis file": "HYS (*.hys)"}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from PySide6.QtWidgets import QStackedWidget, QFileDialog, QMessageBox
    from radia_gui_base import (
        ModePanel, AnalysisWindow, calc_script, msh_output, json_output,
        run_app, COIL_TEMPLATE,
    )
except ImportError as e:
    sys.stderr.write(
        "Radia electromagnet panel requires PySide6 but it could not be "
        "imported:\n"
        "  {}\n\n"
        "Install with:\n"
        "  pip install --upgrade 'radia[gui]'\n".format(e))
    sys.exit(1)


# ============================================================
# Solver display <-> CLI value mappings
# ============================================================
# Consolidated from the old panel's _OMEGA_SOLVERS / _A_SOLVERS /
# _MSC_SOLVERS / _KELVIN_SOLVERS scattered constants.

_FEM_SOLVER_MAP = [
    ("Auto",              "auto"),
    ("Direct (PARDISO)",  "pardiso"),
    ("AMS (Compact)",     "ams"),
    ("BDDC",              "bddc"),
    ("ICCG",              "iccg"),
]
_MSC_SOLVER_MAP = [
    ("LU",       0),
    ("BiCGSTAB", 1),
    ("HACApK",   2),
]


# ============================================================
# Coil-script mixin: "New..." wizard preserved from old panel
# ============================================================

class _CoilTemplateMixin:
    """Adds a coil-template-save wizard button next to the coil_script
    Browse.  Inherited by sub-panels that drive a coil-script-based calc.
    """
    def _on_new_coil_template(self, line_edit):
        existing = line_edit.text().strip()
        if existing and os.path.isdir(os.path.dirname(existing)):
            suggested = os.path.join(os.path.dirname(existing), "coil_new.py")
        else:
            suggested = os.path.abspath("coil_dipole.py")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save new coil template", suggested,
            "Python (*.py);;All (*)")
        if not path:
            return
        try:
            with open(path, "w", encoding="utf-8") as f:
                f.write(COIL_TEMPLATE)
        except OSError as exc:
            QMessageBox.warning(
                self, "Could not write template",
                f"Failed to write {path}:\n{exc}")
            return
        line_edit.setText(path)


# ============================================================
# Sub-panels: one per Formulation, each generator-driven
# ============================================================

class _AccelMagnetPanel(_CoilTemplateMixin, ModePanel):
    """Omega (default) and A-Phi share calc_accel_magnet.py but pass
    --formulation omega/a.  We hide the formulation widget itself since
    the outer Method combo selects it.
    """
    def __init__(self, formulation_cli="omega", parent=None):
        super().__init__(parent)
        from radia.panels.calc_accel_magnet import build_argparser
        self._formulation_cli = formulation_cli
        self.bind_argparser(
            build_argparser(),
            skip=("vol", "output", "formulation", "msh_output"),
            file_browse={
                "coil_script": ("Coil script:", "Python (*.py);;All (*)"),
                "bh_file":     ("BH file:",     "Text files (*.txt *.csv);;All (*)"),
                "hys_file":    ("Hysteresis file:", "HYS files (*.hys);;All (*)"),
            },
            choice_map={"solver": _FEM_SOLVER_MAP},
            labels={
                "fes_order": "FES order:",
                "tol":       "Tolerance:",
                "max_iter":  "Max iterations:",
                "n_steps":   "N steps (FEM):",
                "relax":     "Relaxation:",
            },
        )
        # Preserve the "New..." coil-template wizard button.
        self.add_browse_action(
            "coil_script", "New...", self._on_new_coil_template,
            fixed_width=60)

    def is_runnable(self):
        w = self._widgets.get("coil_script")
        return bool(w and w.text().strip())

    def build_command(self, vol_path):
        coil = self.val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        stem = vol_path or coil
        return self.build_command_from_parser(
            vol_path=vol_path,
            vol_flag="--vol",
            script_path=calc_script("calc_accel_magnet.py"),
            output_path=json_output(stem, "_emfem"),
            extra=["--formulation", self._formulation_cli,
                   "--msh-output", msh_output(stem, "_emfem")],
        )


class _MSCPanel(_CoilTemplateMixin, ModePanel):
    """MSC (Radia hex via calc_accel_msc.py)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        from radia.panels.calc_accel_msc import build_argparser
        self.bind_argparser(
            build_argparser(),
            skip=("vol", "output", "msh_output"),
            file_browse={
                "coil_script": ("Coil script:", "Python (*.py);;All (*)"),
                "bh_file":     ("BH file:",     "Text files (*.txt *.csv);;All (*)"),
                "hys_file":    ("Hysteresis file:", "HYS files (*.hys);;All (*)"),
            },
            choice_map={"solver": _MSC_SOLVER_MAP},
            labels={
                "ima":      "IMA symmetry:",
                "tol":      "Tolerance:",
                "max_iter": "Max iterations:",
                "relax":    "Relaxation:",
            },
        )
        self.add_browse_action(
            "coil_script", "New...", self._on_new_coil_template,
            fixed_width=60)

    def is_runnable(self):
        w = self._widgets.get("coil_script")
        return bool(w and w.text().strip())

    def build_command(self, vol_path):
        coil = self.val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        stem = vol_path or coil
        return self.build_command_from_parser(
            vol_path=vol_path,
            vol_flag="--vol",
            script_path=calc_script("calc_accel_msc.py"),
            output_path=json_output(stem, "_msc"),
            extra=["--msh-output", msh_output(stem, "_msc")],
        )


class _KelvinBenchPanel(ModePanel):
    """Kelvin Benchmark verifies the Cubit-Kelvin pipeline against the
    analytical magnetic-sphere-in-uniform-field solution.  No coil
    script; uniform external field is the source."""
    def __init__(self, parent=None):
        super().__init__(parent)
        from radia.panels.calc_kelvin_benchmark import build_argparser
        self.bind_argparser(
            build_argparser(),
            skip=("vol", "msh_output",
                  "probe_x", "probe_y", "probe_z"),
            labels={
                "fes_order": "FES order:",
                "mu_r":      "mu_r:",
                "H0":        "H0 [A/m]:",
                "field_axis": "Field axis:",
                "R_kelvin":  "Kelvin sphere R [m]:",
            },
        )

    def is_runnable(self):
        return True   # vol_path is required; checked at build_command time

    def build_command(self, vol_path):
        if not vol_path:
            raise ValueError(
                "Kelvin Benchmark requires a .vol mesh with sphere / "
                "kelvin_int / kelvin_ext sidesets.")
        return self.build_command_from_parser(
            vol_path=vol_path,
            vol_flag="--vol",
            script_path=calc_script("calc_kelvin_benchmark.py"),
            output_path=json_output(vol_path, "_kelvin_bench"),
            extra=["--msh-output", msh_output(vol_path, "_kelvin_bench")],
        )


class _ClebschPanel(ModePanel):
    """Clebsch hodograph (reduced-potential) forward mode.  Self-meshes a
    verified 2-D / axisymmetric geometry (cylinder / sphere) in a uniform
    applied field, solves the reduced scalar potential Omega, and
    post-processes the field into its conjugate Clebsch potentials -- the
    flux function psi and the scalar potential Phi -- reporting the
    bidirectional consistency (an a-posteriori field-quality metric) and a
    flux-line / equipotential .msh.  No coil script, no wp_vol: the mode
    self-meshes the verified demonstration geometry.

    The 3-D inverse pole-design (geometry-as-unknown hodograph PDE) is a
    separate research line and is NOT exposed here (repo-first: it stays on
    the branch until verified on a real device)."""
    def __init__(self, parent=None):
        super().__init__(parent)
        from radia.panels.calc_clebsch_hodograph import build_argparser
        self.bind_argparser(
            build_argparser(),
            skip=("output", "msh_output"),     # no "vol" -- self-meshes
            labels={
                "geometry":  "Geometry:",
                "mu_r":      "mu_r:",
                "maxh":      "Mesh size maxh:",
                "fes_order": "FES order:",
            },
        )

    def is_runnable(self):
        return True       # self-meshes; no wp_vol required

    def build_command(self, vol_path):
        # Self-meshing: vol_path is ignored; outputs are keyed on the
        # geometry name (no input file to anchor next to).
        geom = self.val("geometry") or "cylinder"
        return self.build_command_from_parser(
            vol_path="",
            vol_flag=None,
            script_path=calc_script("calc_clebsch_hodograph.py"),
            output_path=json_output(geom, "_clebsch"),
            extra=["--msh-output", msh_output(geom, "_clebsch")],
        )


# ============================================================
# Top-level EMPanel: wp_vol + Method combo + QStackedWidget
# ============================================================

# Method label -> sub-panel factory (eager construction preserves
# per-method widget state across switches).
FORM_OMEGA = "Omega"
FORM_APHI  = "A-Phi"
FORM_MSC   = "MSC"
FORM_KELVIN_BENCH = "Kelvin Benchmark"
FORM_CLEBSCH = "Clebsch hodograph"

_METHOD_FACTORIES = [
    (FORM_OMEGA,        lambda: _AccelMagnetPanel("omega")),
    (FORM_APHI,         lambda: _AccelMagnetPanel("a")),
    (FORM_MSC,          lambda: _MSCPanel()),
    (FORM_KELVIN_BENCH, lambda: _KelvinBenchPanel()),
    (FORM_CLEBSCH,      lambda: _ClebschPanel()),
]


class EMPanel(ModePanel):
    """Composite top-level panel.

    Owns: wp_vol Browse (per-panel .vol, 4.35.0+ legacy mechanism) +
    Method combo + QStackedWidget hosting the 5 sub-panels.

    Per-method state is namespaced as "<method>/<key>" in save_state /
    restore_state, so switching methods preserves what the user typed
    under each.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        # 1. Per-panel mesh .vol Browse (legacy contract for AnalysisWindow.wp_vol_path).
        self.add_browse(
            "wp_vol", "Mesh .vol:",
            filter_str="Netgen Vol (*.vol);;All (*)")

        # 2. Method combo.
        self._method_combo = self.add_combo(
            "method", "Method:",
            [name for name, _ in _METHOD_FACTORIES])

        # 3. QStackedWidget hosting the 5 sub-panels.
        self._sub_panels = {name: factory()
                            for name, factory in _METHOD_FACTORIES}
        self._stack = QStackedWidget()
        for name, _ in _METHOD_FACTORIES:
            self._stack.addWidget(self._sub_panels[name])
        # Add the stack as a single full-width row.
        self._form.addRow(self._stack)

        self._method_combo.currentTextChanged.connect(self._on_method_changed)
        self._on_method_changed(FORM_OMEGA)

    # ---- delegation to current sub-panel ------------------------

    def _current_sub(self):
        return self._sub_panels[self._method_combo.currentText()]

    def _on_method_changed(self, text):
        if text in self._sub_panels:
            self._stack.setCurrentWidget(self._sub_panels[text])
        cb = getattr(self, "validationChanged", None)
        if callable(cb):
            cb()

    def is_runnable(self):
        return self._current_sub().is_runnable()

    def build_command(self, vol_path):
        return self._current_sub().build_command(vol_path)

    def wp_vol_path(self):
        """Source the mesh .vol from the panel's own Browse field (legacy
        4.35.0+ contract honoured by AnalysisWindow._on_run)."""
        return self.val("wp_vol") if "wp_vol" in self._widgets else ""

    # ---- save / restore: per-method namespace -------------------

    def save_state(self):
        # Top-level widgets (wp_vol + method).
        state = super().save_state()
        # Each sub-panel under its method prefix.
        for name, sub in self._sub_panels.items():
            for k, v in sub.save_state().items():
                state[f"{name}/{k}"] = v
        return state

    def restore_state(self, state):
        if not state:
            return
        # Top-level (un-prefixed) keys first.
        top = {k: v for k, v in state.items() if "/" not in k}
        super().restore_state(top)
        # Per-method dispatch.
        for name, sub in self._sub_panels.items():
            prefix = f"{name}/"
            sub.restore_state({k[len(prefix):]: v
                               for k, v in state.items()
                               if k.startswith(prefix)})
        # Update visibility after method combo has been restored.
        self._on_method_changed(self._method_combo.currentText())


# ============================================================
# Window
# ============================================================

class EMWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__("Radia - Electromagnet", vol_path,
                         settings_key="em")
        panel = EMPanel()
        self._set_panel(panel)
        # Restore previous session's widget state FIRST, then let any
        # explicit constructor-supplied vol_path override the restored
        # wp_vol.  The previous order (setText before restore) silently
        # lost the Cubit Solve menu's current .vol path whenever a
        # different .vol was saved last session, because restore would
        # overwrite the just-set wp_vol text.  Mirror radia_ih.py:1452.
        self._restore_settings()
        if vol_path and "wp_vol" in panel._widgets:
            panel._widgets["wp_vol"].setText(self.display_path(vol_path))
        self._update_run_state()


def main():
    run_app(EMWindow)


if __name__ == "__main__":
    main()
