"""
Radia EM (Electromagnet) analysis window.

Formulations (one panel, switched via combo):
  Omega          FEM scalar potential.  Coil-script driven.
  A-Phi          FEM vector potential.  Coil-script driven.
  MSC            Radia magnetic-surface-charge.  Coil-script driven.
  Kelvin Benchmark   Verify Cubit-Kelvin pipeline against analytical
                     magnetic-sphere-in-uniform-field.  Uniform external
                     field input (no coil-script).

Usage:
    python -m radia.radia_em model.vol
    python radia_em.py model.vol
"""

import sys
import os

from PySide6.QtWidgets import QLabel

TITLE = "Electromagnet"
REQUIRED_LABELS = []
OPTIONAL_LABELS = ["iron", "magnetic", "air", "kelvin", "sphere",
                   "kelvin_int", "kelvin_ext"]
OPTIONAL_FILES = {"Coil script": "Python (*.py)", "Hysteresis file": "HYS (*.hys)"}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, msh_output, json_output,
    run_app, _PYTHON,
)


# ---- Formulation labels (combo entries) -----------------------------------
FORM_OMEGA = "Omega"
FORM_APHI = "A-Phi"
FORM_MSC = "MSC"
FORM_KELVIN_BENCH = "Kelvin Benchmark"
FORMULATIONS = [FORM_OMEGA, FORM_APHI, FORM_MSC, FORM_KELVIN_BENCH]


class EMPanel(ModePanel):
    """EM analysis panel.

    Layout (sections):
        Formulation         choose Omega / A-Phi / MSC / Kelvin Benchmark
        Source              coil_script (FEM/MSC) OR uniform H0 (Kelvin)
        Material            mu_r / BH curve / hysteresis file
        Mesh                FES order, Kelvin sphere radius
        Solver              direct / iterative + tolerance
        Iteration           max_iter, relax, n_steps, Picard / Newton
        Symmetry            IMA string (MSC only)

    Sections collapse via `_set_row_visible` when not relevant for the
    current Formulation -- e.g. coil_script is hidden in Kelvin
    Benchmark mode, n_steps is hidden in MSC mode, etc.

    NOTE: FES order >= 2 is supported and Kelvin BC is correctly
    enforced at any order.  See `memory/feedback_kelvin_p2_works.md`
    + `memory/feedback_cubit_kelvin_p_convergence.md`.  Default
    fes_order=1 is kept for the coil-driven Omega path because the
    Periodic-Omega-reduced single-space formulation in
    `calc_accel_magnet.py` shows non-monotonic p-convergence on the
    em_elf_quarter geometry (a separate formulation issue, NOT a
    Kelvin issue).  Kelvin Benchmark mode defaults fes_order=2 since
    its Dirichlet-lift formulation is verified at p>=2.

    Bundled Kelvin Benchmark samples (verified 2026-04-26 at p=2):
        kelvin_benchmark_sphere_1_2.vol   1/2 model, error +1.07%
        kelvin_benchmark_sphere_1_4.vol   1/4 model, error +0.71%
    The 1/8 build script ships
    (`kelvin_benchmark_sphere_1_8_build.py`) but the .vol does NOT
    -- the corner-GND geometry currently produces Hz~0; see
    `memory/feedback_kelvin_1_8_blocker.md`.
    """

    _OMEGA_SOLVERS = ["Direct (PARDISO)", "AMG (Compact)", "BDDC"]
    _A_SOLVERS = ["Direct (PARDISO)", "AMS (Compact)", "BDDC"]
    _MSC_SOLVERS = ["0 (LU)", "1 (BiCGSTAB)", "2 (HACApK)"]
    _KELVIN_SOLVERS = ["Direct (PARDISO)"]
    _SOLVER_MAP = {
        "Direct (PARDISO)": "pardiso",
        "AMG (Compact)": "amg",
        "AMS (Compact)": "ams",
        "BDDC": "bddc",
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    # ------------------------------------------------------------------
    # Section header helper (copied from radia_ih.py 2026-04-25 -- the
    # section style with explicit fixed-height label and padding is the
    # one that survives both offscreen + real-desktop Qt rendering).
    # ------------------------------------------------------------------
    def _add_section(self, title, key=None):
        """Insert a bold section header row into the form.

        ``key`` registers the row in ``self._row_indices`` so it can be
        collapsed together with the widgets it groups via
        ``_set_row_visible(key, False)``.
        """
        lbl = QLabel(f"<b>{title}</b>")
        f = lbl.font(); f.setPointSize(f.pointSize() + 1); lbl.setFont(f)
        lbl.setStyleSheet("QLabel { color: #444; }")
        fm_h = lbl.fontMetrics().boundingRect("Mg").height()
        lbl.setFixedHeight(fm_h + 10)
        self._form.addRow("", lbl)
        if key is not None:
            self._row_indices[key] = self._form.rowCount() - 1

    # ------------------------------------------------------------------
    # UI build
    # ------------------------------------------------------------------
    def _build_ui(self):
        # ---- Formulation ----
        self._add_section("Formulation")
        self._form_combo = self.add_combo(
            "formulation", "Formulation:", FORMULATIONS)
        self._form_combo.currentTextChanged.connect(
            self._on_formulation_changed)

        # ---- Source ----
        self._add_section("Source", key="_sec_source")
        # FEM / MSC use a coil script.
        self.add_browse("coil_script", "Coil script:",
                        filter_str="Python (*.py);;All (*)")
        # Kelvin Benchmark uses a uniform external field instead.
        self.add_line("H0", "H0 [A/m]:", "1.0")
        self.add_combo("field_axis", "Field axis:", ["x", "y", "z"])
        # Default field axis = z
        self._widgets["field_axis"].setCurrentIndex(2)

        # ---- Material ----
        self._add_section("Material", key="_sec_material")
        mat_combo = self.add_combo(
            "material", "Material:",
            ["mu_r (Linear)", "BH Curve", "Hysteresis (.hys)"])
        self.add_line("mu_r", "mu_r:", "1000")
        self.add_browse("bh_file", "BH file:",
                        filter_str="Text files (*.txt *.csv);;All (*)")
        self.add_browse("hys_file", "Hysteresis file:",
                        filter_str="HYS files (*.hys);;All (*)")
        mat_combo.currentIndexChanged.connect(self._on_material_changed)

        # ---- Mesh ----
        self._add_section("Mesh", key="_sec_mesh")
        self.add_spin("fes_order", "FES order:", 1, 1, 5)
        # Kelvin Benchmark only -- declared here so panel layout stays
        # consistent across modes.  Hidden when not Kelvin Benchmark.
        self.add_line("R_kelvin", "Kelvin sphere R [m]:", "0.20")

        # ---- Solver ----
        self._add_section("Solver", key="_sec_solver")
        self.add_combo("solver", "Solver:", self._OMEGA_SOLVERS)
        self.add_line("tol", "Tolerance:", "1e-3")

        # ---- Iteration ----
        self._add_section("Iteration", key="_sec_iter")
        self.add_spin("max_iter", "Max iterations:", 100, 1, 9999)
        self.add_line("relax", "Relaxation:", "0.0")
        # FEM-specific: nonlinear time stepping.
        self.add_spin("n_steps", "N steps (FEM):", 1, 1, 100)
        self.add_combo("method", "Method (FEM):", ["Picard", "Newton"])

        # ---- Symmetry (MSC only) ----
        self._add_section("Symmetry (MSC)", key="_sec_sym")
        self.add_line("ima", "IMA symmetry:", "",
                      placeholder="auto from sym_*_* (blank) / e.g. +x-z")

        # Initial state: pick Omega defaults.
        self._on_formulation_changed(FORM_OMEGA)
        self._on_material_changed(0)

    # ------------------------------------------------------------------
    # Mode visibility logic
    # ------------------------------------------------------------------
    def _on_formulation_changed(self, text):
        is_msc = (text == FORM_MSC)
        is_kelvin = (text == FORM_KELVIN_BENCH)
        is_fem = not (is_msc or is_kelvin)

        # Source row visibility:
        #   coil_script  -- FEM, MSC (NOT Kelvin Benchmark)
        #   H0, axis     -- Kelvin Benchmark only
        self._set_row_visible("coil_script", is_fem or is_msc)
        self._set_row_visible("H0", is_kelvin)
        self._set_row_visible("field_axis", is_kelvin)

        # Mesh: R_kelvin is Kelvin Benchmark only.
        self._set_row_visible("R_kelvin", is_kelvin)

        # FEM-only knobs.
        self._set_row_visible("fes_order", is_fem or is_kelvin)
        self._set_row_visible("n_steps", is_fem)
        self._set_row_visible("method", is_fem)

        # Iteration: tol/max_iter/relax used by FEM nonlinear and MSC.
        # Kelvin Benchmark is a single linear solve -- hide them.
        self._set_row_visible("tol", not is_kelvin)
        self._set_row_visible("max_iter", not is_kelvin)
        self._set_row_visible("relax", not is_kelvin)

        # Symmetry section: MSC only.
        self._set_row_visible("ima", is_msc)
        self._set_row_visible("_sec_sym", is_msc)

        # Mesh section: hide for MSC (no fes_order, no Kelvin radius).
        self._set_row_visible("_sec_mesh", not is_msc)

        # Iteration section: hide for Kelvin Benchmark (single linear solve).
        self._set_row_visible("_sec_iter", not is_kelvin)

        # Solver list per formulation.
        combo = self._widgets["solver"]
        prev = combo.currentText()
        combo.clear()
        if is_msc:
            combo.addItems(self._MSC_SOLVERS)
        elif text == FORM_APHI:
            combo.addItems(self._A_SOLVERS)
        elif is_kelvin:
            combo.addItems(self._KELVIN_SOLVERS)
        else:
            combo.addItems(self._OMEGA_SOLVERS)
        idx = combo.findText(prev)
        if idx >= 0:
            combo.setCurrentIndex(idx)

        # Default fes_order: 2 for Kelvin Benchmark (verified at p>=2).
        if is_kelvin:
            self._widgets["fes_order"].setValue(2)

        cb = getattr(self, 'validationChanged', None)
        if callable(cb):
            cb()

    def _on_material_changed(self, idx):
        self._set_row_visible("mu_r", idx == 0)
        self._set_row_visible("bh_file", idx == 1)
        self._set_row_visible("hys_file", idx == 2)

    # ------------------------------------------------------------------
    # Run gating
    # ------------------------------------------------------------------
    def is_runnable(self):
        formulation = self.val("formulation")
        if formulation == FORM_KELVIN_BENCH:
            # Kelvin Benchmark needs only a .vol with sphere/kelvin_int/
            # kelvin_ext sidesets.  No coil script.
            return True
        return bool(self._widgets["coil_script"].text().strip())

    # ------------------------------------------------------------------
    # Command builders (one per formulation)
    # ------------------------------------------------------------------
    def build_command(self, vol_path):
        formulation = self.val("formulation")
        if formulation == FORM_MSC:
            return self._build_msc_command(vol_path)
        if formulation == FORM_KELVIN_BENCH:
            return self._build_kelvin_command(vol_path)
        return self._build_fem_command(vol_path, formulation)

    def _material_args(self, mat_keyword_for_linear="linear"):
        """Shared --material/--mu-r/--bh-file/--hys-file argument block."""
        idx = self._widgets["material"].currentIndex()
        if idx == 0:
            return ["--material", mat_keyword_for_linear,
                    "--mu-r", self.val("mu_r")]
        if idx == 1:
            args = ["--material", "steel"]
            bh = self.val("bh_file")
            if bh:
                args += ["--bh-file", bh]
            return args
        args = ["--material", "hysteresis"]
        hys = self.val("hys_file")
        if hys:
            args += ["--hys-file", hys]
        return args

    # CLI-DIFF: ignore --output --sigma -- output auto-generated by panel; sigma defaulted from --material preset
    def _build_fem_command(self, vol_path, formulation):
        coil = self.val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        solver_key = self._SOLVER_MAP.get(self.val("solver"), "pardiso")

        cmd = [_PYTHON, calc_script("calc_accel_magnet.py"),
               "--formulation", formulation.lower(),
               "--fes-order", self.val("fes_order"),
               "--coil-script", coil,
               "--solver", solver_key,
               "--tol", self.val("tol"),
               "--max-iter", self.val("max_iter"),
               "--n-steps", self.val("n_steps"),
               "--relax", self.val("relax")]
        cmd += self._material_args("linear")
        if vol_path:
            cmd += ["--vol", vol_path]
        # Note: --ima is MSC-only (Radia image-method symmetry).  The
        # FEM formulation in calc_accel_magnet.py uses standard FEM
        # BC (sym_tangential / sym_normal sidesets) instead and does
        # not declare --ima, so do not emit it here.
        if self.val("method") == "Newton":
            cmd += ["--newton"]
        cmd += ["--msh-output", msh_output(vol_path or coil, "_emfem"),
                "--output", json_output(vol_path or coil, "_emfem")]
        return cmd

    def _build_msc_command(self, vol_path):
        coil = self.val("coil_script")
        if not coil:
            raise ValueError("No coil script specified.")
        solver_id = self.val("solver").split()[0]

        cmd = [_PYTHON, calc_script("calc_accel_msc.py"),
               "--coil-script", coil,
               "--solver", solver_id,
               "--tol", self.val("tol"),
               "--max-iter", self.val("max_iter"),
               "--relax", self.val("relax")]
        cmd += self._material_args("custom")
        if vol_path:
            cmd += ["--vol", vol_path]
        ima = self.val("ima")
        if ima:
            cmd += ["--ima", ima]
        cmd += ["--msh-output", msh_output(vol_path or coil, "_msc"),
                "--output", json_output(vol_path or coil, "_msc")]
        return cmd

    # CLI-DIFF: ignore --probe-x --probe-y --probe-z -- panel always probes the origin (verification default)
    def _build_kelvin_command(self, vol_path):
        if not vol_path:
            raise ValueError(
                "Kelvin Benchmark requires a .vol mesh with sphere / "
                "kelvin_int / kelvin_ext sidesets.")
        # Kelvin Benchmark always uses the linear material (mu_r),
        # since the analytical reference assumes a linear sphere.
        mu_r = self.val("mu_r") or "100"
        cmd = [_PYTHON, calc_script("calc_kelvin_benchmark.py"),
               "--vol", vol_path,
               "--fes-order", self.val("fes_order"),
               "--mu-r", mu_r,
               "--H0", self.val("H0"),
               "--field-axis", self.val("field_axis"),
               "--R-kelvin", self.val("R_kelvin")]
        cmd += ["--msh-output", msh_output(vol_path, "_kelvin_bench"),
                "--output", json_output(vol_path, "_kelvin_bench")]
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
