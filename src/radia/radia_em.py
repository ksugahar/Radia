"""
Radia EM (Electromagnet) analysis window.

Modes: A-Phi (FEM) / Omega (FEM) / MSC (Radia)
Switch via combo box -- single window.

Usage:
    python -m radia.radia_em model.vol
    python radia_em.py model.vol
"""

import sys
import os

TITLE = "Electromagnet"
REQUIRED_LABELS = []
OPTIONAL_LABELS = ["iron", "air", "kelvin_int", "kelvin_ext"]
OPTIONAL_FILES = {"Coil script": "Python (*.py)", "Hysteresis file": "HYS (*.hys)"}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, msh_output, run_app, _PYTHON,
)


class EMPanel(ModePanel):

    _OMEGA_SOLVERS = ["Direct (PARDISO)", "AMG (Compact)", "BDDC"]
    _A_SOLVERS = ["Direct (PARDISO)", "AMS (Compact)", "BDDC"]
    _MSC_SOLVERS = ["0 (LU)", "1 (BiCGSTAB)", "2 (HACApK)"]
    _SOLVER_MAP = {
        "Direct (PARDISO)": "pardiso",
        "AMG (Compact)": "amg",
        "AMS (Compact)": "ams",
        "BDDC": "bddc",
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        # Formulation selector
        self._form_combo = self.add_combo(
            "formulation", "Formulation:", ["Omega", "A-Phi", "MSC"])
        self._form_combo.currentTextChanged.connect(
            self._on_formulation_changed)

        # FEM fields
        self.add_spin("fes_order", "FES order:", 1, 1, 5)
        self.add_browse("coil_script", "Coil script:",
                        filter_str="Python (*.py);;All (*)")

        # Material: mu_r / BH / Hysteresis
        mat_combo = self.add_combo("material", "Material:",
                                   ["mu_r (Linear)", "BH Curve",
                                    "Hysteresis (.hys)"])
        self.add_line("mu_r", "mu_r:", "1000")
        self.add_browse("bh_file", "BH file:",
                        filter_str="Text files (*.txt *.csv);;All (*)")
        self.add_browse("hys_file", "Hysteresis file:",
                        filter_str="HYS files (*.hys);;All (*)")
        mat_combo.currentIndexChanged.connect(self._on_material_changed)

        # Solver (content depends on formulation)
        self.add_combo("solver", "Solver:", self._OMEGA_SOLVERS)

        # MSC-only: Radia image-method symmetry string.  Leave blank
        # to let calc_accel_msc auto-detect from `sym_<bc>_<axis>`
        # sidesets in the .vol (sym_bn=0_x -> '+x', sym_ht=0_x ->
        # '-x', combined in x/y/z order).  Pass an explicit string to
        # override.  (FEM formulation picks Dirichlet/Natural per
        # formulation directly from the same sidesets; no widget.)
        self.add_line("ima", "IMA symmetry:", "",
                      placeholder="auto from sym_*_* (blank) / e.g. +x-z")
        self.add_line("tol", "Tolerance:", "1e-3")
        self.add_spin("max_iter", "Max iterations:", 100, 1, 9999)
        self.add_line("relax", "Relaxation:", "0.0")

        # FEM-only
        self.add_spin("n_steps", "N steps:", 1, 1, 100)
        self.add_combo("method", "Method:", ["Picard", "Newton"])

        # (--unit-scale widget was retired 2026-04-24 with the
        # cub5 -> vol migration; radia_export netgen always writes
        # coordinates in meters.  Panel keeps no widget for it.)

        # Initial state
        self._on_formulation_changed("Omega")
        self._on_material_changed(0)

    def _on_formulation_changed(self, text):
        is_msc = (text == "MSC")
        is_fem = not is_msc

        # FEM-only fields
        for key in ("fes_order", "n_steps", "method"):
            self._set_row_visible(key, is_fem)

        # (unit_scale widget removed along with --unit-scale CLI)
        # IMA is MSC-only (see widget comment above).
        self._set_row_visible("ima", is_msc)

        # Solver list
        combo = self._widgets["solver"]
        prev = combo.currentText()
        combo.clear()
        if is_msc:
            combo.addItems(self._MSC_SOLVERS)
        elif text == "A-Phi":
            combo.addItems(self._A_SOLVERS)
        else:
            combo.addItems(self._OMEGA_SOLVERS)
        idx = combo.findText(prev)
        if idx >= 0:
            combo.setCurrentIndex(idx)

        cb = getattr(self, 'validationChanged', None)
        if callable(cb):
            cb()

    def _on_material_changed(self, idx):
        self._set_row_visible("mu_r", idx == 0)
        self._set_row_visible("bh_file", idx == 1)
        self._set_row_visible("hys_file", idx == 2)

    def is_runnable(self):
        return bool(self._widgets["coil_script"].text().strip())

    def build_command(self, vol_path):
        coil = self.val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        formulation = self.val("formulation")

        if formulation == "MSC":
            return self._build_msc_command(vol_path, coil)
        else:
            return self._build_fem_command(vol_path, coil, formulation)

    # CLI-DIFF: ignore --output --sigma -- output auto-generated by panel; sigma defaulted from --material preset
    def _build_fem_command(self, vol_path, coil, formulation):
        solver_key = self._SOLVER_MAP.get(self.val("solver"), "pardiso")

        mat_idx = self._widgets["material"].currentIndex()
        if mat_idx == 0:
            mat_args = ["--material", "linear",
                        "--mu-r", self.val("mu_r")]
        elif mat_idx == 1:
            mat_args = ["--material", "steel"]
            bh = self.val("bh_file")
            if bh:
                mat_args += ["--bh-file", bh]
        else:
            mat_args = ["--material", "hysteresis"]
            hys = self.val("hys_file")
            if hys:
                mat_args += ["--hys-file", hys]

        cmd = [_PYTHON, calc_script("calc_accel_magnet.py"),
               "--formulation", formulation.lower(),
               "--fes-order", self.val("fes_order"),
               "--coil-script", coil,
               "--solver", solver_key,
               "--tol", self.val("tol"),
               "--max-iter", self.val("max_iter"),
               "--n-steps", self.val("n_steps"),
               "--relax", self.val("relax")]
        cmd += mat_args
        if vol_path:
            cmd += ["--vol", vol_path]
        # Note: --ima is MSC-only (Radia image-method symmetry).  The
        # FEM formulation in calc_accel_magnet.py uses standard FEM
        # BC (sym_tangential / sym_normal sidesets) instead and does
        # not declare --ima, so do not emit it here.
        if self.val("method") == "Newton":
            cmd += ["--newton"]
        cmd += ["--msh-output", msh_output(
            vol_path or coil, "_emfem")]
        return cmd

    def _build_msc_command(self, vol_path, coil):
        solver_id = self.val("solver").split()[0]

        mat_idx = self._widgets["material"].currentIndex()
        if mat_idx == 0:
            mat_args = ["--material", "custom",
                        "--mu-r", self.val("mu_r")]
        elif mat_idx == 1:
            mat_args = ["--material", "steel"]
            bh = self.val("bh_file")
            if bh:
                mat_args += ["--bh-file", bh]
        else:
            mat_args = ["--material", "hysteresis"]
            hys = self.val("hys_file")
            if hys:
                mat_args += ["--hys-file", hys]

        cmd = [_PYTHON, calc_script("calc_accel_msc.py"),
               "--coil-script", coil,
               "--solver", solver_id,
               "--tol", self.val("tol"),
               "--max-iter", self.val("max_iter"),
               "--relax", self.val("relax")]
        # Note: --unit-scale was retired 2026-04-24 when the panel
        # migrated from .cub5 to .vol (radia_export netgen always
        # writes coordinates in meters).
        cmd += mat_args
        if vol_path:
            cmd += ["--vol", vol_path]
        ima = self.val("ima")
        if ima:
            cmd += ["--ima", ima]
        cmd += ["--msh-output", msh_output(
            vol_path or coil, "_msc")]
        return cmd


class EMWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__("Radia - Electromagnet", vol_path,
                         settings_key="em")
        panel = EMPanel()
        self._set_panel(panel)
        self._restore_settings()
        self._update_run_state()


def main():
    run_app(EMWindow)


if __name__ == "__main__":
    main()
