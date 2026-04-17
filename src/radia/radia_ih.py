"""
Radia IH (Induction Heating) analysis window.

Modes:
  PEEC     -- Filament coil (no mesh needed for coil), PRIMA broadband.
              Fastest for coil self-impedance / freq sweep.
  PEEC+FEM -- PEEC coil + BEM-SIBC workpiece coupling (Picard, 1 step).
              .vol needed only for workpiece surface mesh.
  BEM      -- EFIE source/sink coil + optional workpiece SIBC / ESIM.
              Legacy; slower than PEEC for coil inductance.
  FEM      -- Kelvin + SIBC / ESIM (full volume mesh).

Switch via combo box -- single window.

Usage:
    python -m radia.radia_ih model.vol
    python radia_ih.py model.vol
"""

import sys
import os

TITLE = "Induction Heating"
# coil:   block name (volume label) — used by calc_inductance.py to filter
#         the BEM surface mesh to coil-adjacent faces only.
# source: sideset name on one coil terminal face (current injection).
# sink:   sideset name on the other coil terminal face (current extraction).
# BEM requires source/sink sidesets; PEEC does not (coil defined parametrically).
# The label check is mode-dependent; see is_runnable().
REQUIRED_LABELS = []   # relaxed — PEEC needs no mesh labels at all
OPTIONAL_LABELS = ["workpiece", "air", "kelvin", "kelvin_int", "kelvin_ext"]
OPTIONAL_FILES = {}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, msh_output, run_app, _PYTHON,
)


class IHPanel(ModePanel):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._build_ui()

    def _build_ui(self):
        # Method selector — defaults to BEM. user_settings restored
        # later (in IHWindow.__init__) overrides this when present.
        self._method_combo = self.add_combo(
            "method", "Method:", ["PEEC+FEM", "FEM"])
        self._method_combo.currentTextChanged.connect(self._on_method_changed)

        self.add_spin("fes_order", "FES order:", 0, 0, 5)

        # Solver row — context-dependent items, populated by
        # _on_method_changed. Single combo so the layout does NOT
        # shift when switching modes.
        solver_combo = self.add_combo("solver", "Solver:", ["LU (direct)"])

        # FEM-only iteration cap. Hidden in BEM mode (BEM is direct).
        self.add_spin("max_iter", "Max iterations:", 15, 1, 200)

        # ============ Coil parameters ============
        # Frequency, current, sigma — applicable to both BEM and FEM.
        self.add_line("freq", "Frequency [Hz]:", "50000")
        self.add_line("current", "Coil current [A]:", "1.0")
        self.add_line("coil_sigma", "Coil sigma [S/m]:", "5.8e7")

        # ============ PEEC coil parameters (no mesh needed) ============
        self.add_line("coil_radius", "Coil radius [m]:", "0.05")
        self.add_line("coil_pitch", "Coil pitch [m]:", "0.008")
        self.add_line("coil_turns", "Turns:", "5")
        self.add_line("wire_w", "Wire width [m]:", "0.006")
        self.add_line("wire_h", "Wire height [m]:", "0.006")
        self.add_spin("nwinc", "nwinc (skin subdiv):", 3, 1, 10)
        self.add_spin("seg_per_turn", "Segments/turn:", 20, 8, 50)
        self.add_spin("prima_q", "PRIMA order:", 30, 5, 100)

        # ============ Workpiece parameters ============
        # The same combo drives both BEM and FEM:
        #   off  -> workpiece NOT in DOF (coil-only inductance, BEM only)
        #   SIBC -> linear surface impedance + per-panel curvature
        #   ESIM -> nonlinear 1D cell problem + per-panel curvature
        wp_combo = self.add_combo("workpiece_mode", "Workpiece:",
                                  ["off", "SIBC", "ESIM"])
        wp_combo.currentTextChanged.connect(self._on_workpiece_changed)

        self.add_line("wp_sigma", "WP sigma [S/m]:", "2e6")
        self.add_line("half_thickness", "Half thickness [m]:", "0.005")
        self.add_line("mu_r", "mu_r:", "100")
        self.add_browse("bh_file", "BH file:",
                        filter_str="Text files (*.txt *.csv);;All (*)")
        # ESIM 1D cell-problem coordinate system.
        self.add_combo("esim_geometry", "ESIM geometry:",
                       ["local_curvature", "none"])

        # Air field calc — BEM only (FEM volume mesh always solves
        # for the field everywhere, so the toggle is meaningless in
        # FEM and would mislead the user). Hidden in FEM mode.
        self.add_combo("air_mode", "Air field calc:", ["off", "on"])

        # Initial state — default Method = PEEC.
        self._method_combo.setCurrentText("PEEC+FEM")
        self._on_method_changed("PEEC+FEM")
        self._on_validation_changed()

    def _on_method_changed(self, method):
        is_peec = (method == "PEEC+FEM")
        is_fem = (method == "FEM")

        # Solver combo: same widget, different items per method.
        solver = self._widgets["solver"]
        prev = solver.currentText()
        solver.clear()
        if is_peec:
            solver.addItems(["Dense LU", "HACApK saddle", "PRIMA"])
        else:
            solver.addItems(["pardiso", "bddc", "iccg", "ams"])
        idx = solver.findText(prev)
        if idx >= 0:
            solver.setCurrentIndex(idx)

        # PEEC coil parameters — visible only in PEEC modes.
        for key in ("coil_radius", "coil_pitch", "coil_turns",
                     "wire_w", "wire_h", "nwinc", "seg_per_turn"):
            self._set_row_visible(key, is_peec)
        self._set_row_visible("prima_q",
                              is_peec and self.val("solver") == "PRIMA")

        # FEM-only widgets.
        self._set_row_visible("max_iter", is_fem)
        self._set_row_visible("fes_order", is_fem)
        self._set_row_visible("air_mode", False)  # BEM removed

        # Workpiece: always visible. "off" = coil only; SIBC/ESIM = coupling.
        # FEM requires workpiece (forces non-off).
        self._set_row_visible("workpiece_mode", True)
        if is_fem and self.val("workpiece_mode") == "off":
            wp = self._widgets["workpiece_mode"]
            sibc_idx = wp.findText("SIBC")
            if sibc_idx >= 0:
                wp.setCurrentIndex(sibc_idx)

        self._on_workpiece_changed()
        self._on_validation_changed()

    def _on_workpiece_changed(self, _text=None):
        wp_mode = self.val("workpiece_mode")
        has_wp = (wp_mode != "off")

        # Hide all WP-detail widgets first, then show only the
        # subset that matches the current SIBC vs ESIM choice.
        for key in ("wp_sigma", "half_thickness", "mu_r",
                     "bh_file", "esim_geometry"):
            self._set_row_visible(key, False)

        if has_wp:
            is_esim = (wp_mode == "ESIM")
            self._set_row_visible("wp_sigma", True)
            self._set_row_visible("half_thickness", True)
            self._set_row_visible("mu_r", not is_esim)
            self._set_row_visible("bh_file", is_esim)
            self._set_row_visible("esim_geometry", is_esim)

    def _on_validation_changed(self, _text=None):
        cb = getattr(self, 'validationChanged', None)
        if callable(cb):
            cb()

    def is_runnable(self):
        # All methods rely on .jou naming convention; no GUI label inputs.
        return True

    def build_command(self, vol_path):
        method = self.val("method")
        if method == "PEEC+FEM":
            return self._build_peec_command(vol_path)
        if not vol_path:
            raise ValueError("No .vol file specified.")
        return self._build_fem_command(vol_path)

    # PEEC solver combo -> calc-script --solver value
    _PEEC_SOLVER_MAP = {
        "PRIMA": "prima",
        "HACApK saddle": "hacapk",
        "Dense LU": "dense",
    }

    def _build_peec_command(self, vol_path):
        cmd = [_PYTHON, calc_script("calc_peec.py"),
               "--coil-radius", self.val("coil_radius"),
               "--coil-pitch", self.val("coil_pitch"),
               "--coil-turns", self.val("coil_turns"),
               "--wire-w", self.val("wire_w"),
               "--wire-h", self.val("wire_h"),
               "--nwinc", str(self.val("nwinc")),
               "--seg-per-turn", str(self.val("seg_per_turn")),
               "--frequency", self.val("freq"),
               "--current", self.val("current"),
               "--coil-sigma", self.val("coil_sigma"),
               "--solver",
               self._PEEC_SOLVER_MAP.get(self.val("solver"), "dense"),
               "--prima-q", str(self.val("prima_q"))]
        # Workpiece coupling: "off" = coil only, SIBC/ESIM = BEM coupling
        wp_mode = self.val("workpiece_mode")
        if wp_mode != "off" and vol_path:
            imp = "esim" if wp_mode == "ESIM" else "sibc"
            cmd += ["--vol", vol_path,
                    "--workpiece", "sibc",
                    "--impedance-model", imp,
                    "--sigma", self.val("wp_sigma"),
                    "--half-thickness", self.val("half_thickness"),
                    "--mu-r", self.val("mu_r"),
                    "--use-local-curvature"]
        return cmd

    # Solver combo display name -> calc-script --solver value
    _BEM_SOLVER_MAP = {
        "LU (direct)": "lu",
        "MINRES": "minres",
        "GMRES": "gmres",
    }

    def _build_bem_command(self, vol_path):
        # Source/sink/coil labels follow the .jou naming convention
        # ("source"/"sink"/"coil") and are picked up by
        # calc_inductance.py defaults. No label arguments from the GUI.
        cmd = [_PYTHON, calc_script("calc_inductance.py"),
               "--vol", vol_path,
               "--frequency", self.val("freq"),
               "--current", self.val("current"),
               "--msh-output", msh_output(vol_path, "_bem")]
        coil_sigma = self.val("coil_sigma")
        if coil_sigma:
            cmd += ["--coil-sigma", coil_sigma]
        fes = self.val("fes_order")
        if fes and fes != "0":
            cmd += ["--fes-order", fes]

        cmd += ["--solver",
                self._BEM_SOLVER_MAP.get(self.val("solver"), "lu")]

        # Workpiece coupling: combo -> --workpiece + impedance.
        # The workpiece sideset is named "sibc" by .jou convention
        # (see ih_bem_sample.jou).
        wp_mode = self.val("workpiece_mode")
        if wp_mode != "off":
            imp = "esim" if wp_mode == "ESIM" else "sibc"
            cmd += ["--workpiece", "sibc",
                    "--impedance-model", imp,
                    "--sigma", self.val("wp_sigma"),
                    "--half-thickness", self.val("half_thickness"),
                    "--mu-r", self.val("mu_r")]
            if imp == "esim":
                cmd += ["--esim-geometry", self.val("esim_geometry")]
                bh = self.val("bh_file")
                if bh:
                    cmd += ["--bh-file", bh]
            # Per-panel local curvature is now always on (the global
            # half_thickness is only used as a fallback when the
            # extractor cannot recover the local radius).
            cmd += ["--use-local-curvature"]

        # Air field post-processing on/off (BEM only)
        if self.val("air_mode") == "on":
            cmd += ["--field-air"]

        return cmd

    def _build_fem_command(self, vol_path):
        # The same widget set drives FEM. wp_sigma + mu_r come from
        # the shared workpiece group, so the FEM command no longer
        # needs the legacy linear / BH-curve material combo.
        wp_mode = self.val("workpiece_mode")
        impedance = "esim" if wp_mode == "ESIM" else "sibc"
        cmd = [_PYTHON, calc_script("calc_fem_kelvin.py"),
               "--vol", vol_path,
               "--fes-order", self.val("fes_order"),
               "--frequency", self.val("freq"),
               "--current", self.val("current"),
               "--material", "custom",
               "--sigma", self.val("wp_sigma"),
               "--mu-r", self.val("mu_r"),
               "--half-thickness", self.val("half_thickness"),
               "--impedance", impedance,
               "--solver", self.val("solver"),
               "--max-iter", self.val("max_iter"),
               "--msh-output", msh_output(vol_path, "_fem")]
        if impedance == "esim":
            bh = self.val("bh_file")
            if bh:
                cmd += ["--bh-file", bh]
        return cmd


class IHWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__("Radia - Induction Heating", vol_path,
                         settings_key="ih")
        panel = IHPanel()
        self._set_panel(panel)
        self._restore_settings()
        # Auto-populate workpiece/air fields from the .vol's materials
        # so the user does not have to retype labels that the .jou
        # already declared.
        self._populate_optional_labels(vol_path)
        self._update_run_state()

    def _populate_optional_labels(self, vol_path):
        """Read materials from the .vol and pre-set workpiece/air combos
        when the corresponding labels exist. The user can still flip the
        combo back to "off" to disable that feature.
        """
        if not vol_path or not os.path.isfile(vol_path):
            return
        from netgen.meshing import Mesh as NgMesh
        ng = NgMesh()
        ng.Load(vol_path)
        mats = set()
        for i in range(1, 64):
            try:
                n = ng.GetMaterial(i)
            except Exception:
                break
            if n:
                mats.add(n)
        widgets = self._panel._widgets
        # workpiece combo: default to SIBC when a "workpiece" material
        # is present in the .vol (the user can flip to ESIM or off).
        if "workpiece" in mats and "workpiece_mode" in widgets:
            wp = widgets["workpiece_mode"]
            if wp.currentText() == "off":
                idx = wp.findText("SIBC")
                if idx >= 0:
                    wp.setCurrentIndex(idx)
        # air combo: default to "on" when an "air" material exists.
        if "air" in mats and "air_mode" in widgets:
            ac = widgets["air_mode"]
            if ac.currentText() == "off":
                idx = ac.findText("on")
                if idx >= 0:
                    ac.setCurrentIndex(idx)


def main():
    run_app(IHWindow)


if __name__ == "__main__":
    main()
