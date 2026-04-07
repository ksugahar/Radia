"""
Radia IH (Induction Heating) analysis window.

Modes: FEM (Kelvin + SIBC/ESIM) / BEM (EFIE + SIBC)
Switch via combo box -- single window.

Usage:
    python -m radia.radia_ih model.vol
    python radia_ih.py model.vol
"""

import sys
import os

TITLE = "Induction Heating"
REQUIRED_LABELS = ["source", "sink"]
OPTIONAL_LABELS = ["workpiece", "air"]
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
        # Method selector
        self._method_combo = self.add_combo(
            "method", "Method:", ["BEM", "FEM"])
        self._method_combo.currentTextChanged.connect(self._on_method_changed)

        self.add_spin("fes_order", "FES order:", 0, 0, 5)

        # Source / Sink (required for BEM, red if empty)
        src_w = self.add_line("source", "Source block:", "",
                              placeholder="e.g. source")
        sink_w = self.add_line("sink", "Sink block:", "",
                               placeholder="e.g. sink")
        src_w.textChanged.connect(self._on_validation_changed)
        sink_w.textChanged.connect(self._on_validation_changed)

        # Coil parameters (always visible)
        self.add_line("freq", "Frequency [Hz]:", "50000")
        self.add_line("coil_sigma", "Coil sigma [S/m]:", "5.8e7")

        # Workpiece (optional for BEM, always for FEM)
        wp_w = self.add_line("workpiece", "Workpiece block:", "",
                             placeholder="(empty = inductance only)")
        wp_w.textChanged.connect(self._on_workpiece_changed)

        # Workpiece impedance model
        imp_combo = self.add_combo("impedance", "Impedance model:",
                                   ["dowell", "esim", "bem-sibc", "sibc"])
        imp_combo.currentTextChanged.connect(self._on_impedance_changed)

        self.add_line("wp_sigma", "WP sigma [S/m]:", "2e6")
        self.add_line("half_thickness", "Half thickness [m]:", "0.005")
        self.add_line("mu_r", "mu_r:", "100")
        self.add_browse("bh_file", "BH file:",
                        filter_str="Text files (*.txt *.csv);;All (*)")
        self.add_combo("esim_geometry", "ESIM geometry:",
                       ["local_curvature", "planar"])

        # FEM-specific
        self.add_line("current", "Current [A]:", "1.0")
        self.add_line("a_coil", "Coil radius [m]:", "0.003")
        self.add_line("r_wp", "Workpiece radius [m]:", "0.010")
        self.add_combo("solver", "Solver:",
                       ["pardiso", "bddc", "iccg", "ams"])
        self.add_spin("max_iter", "Max iterations:", 15, 1, 200)

        # BEM-specific: air domain
        self.add_line("air", "Air block:", "",
                      placeholder="(empty = no field calculation)")

        # Material for FEM
        mat_combo = self.add_combo("fem_material", "Material:",
                                   ["mu_r (Linear)", "BH Curve"])
        mat_combo.currentIndexChanged.connect(self._on_fem_material_changed)
        self.add_browse("fem_bh_file", "BH file (FEM):",
                        filter_str="Text files (*.txt *.csv);;All (*)")

        # Initial state
        self._on_method_changed("BEM")
        self._on_validation_changed()

    def _on_method_changed(self, method):
        is_bem = (method == "BEM")
        is_fem = not is_bem

        # BEM fields
        for key in ("source", "sink", "air"):
            self._set_row_visible(key, is_bem)

        # FEM fields
        for key in ("current", "a_coil", "r_wp", "solver",
                     "max_iter", "fem_material", "fem_bh_file"):
            self._set_row_visible(key, is_fem)

        # Impedance options differ
        combo = self._widgets["impedance"]
        prev = combo.currentText()
        combo.clear()
        if is_bem:
            combo.addItems(["dowell", "esim", "bem-sibc"])
        else:
            combo.addItems(["sibc", "esim"])
        # Try to keep selection
        idx = combo.findText(prev)
        if idx >= 0:
            combo.setCurrentIndex(idx)

        self._on_workpiece_changed()
        self._on_validation_changed()

    def _on_workpiece_changed(self, _text=None):
        has_wp = bool(self._widgets["workpiece"].text().strip())
        method = self.val("method")

        # Hide all wp params first
        for key in ("impedance", "wp_sigma", "half_thickness",
                     "mu_r", "bh_file", "esim_geometry"):
            self._set_row_visible(key, False)

        if has_wp:
            self._set_row_visible("impedance", True)
            self._set_row_visible("wp_sigma", True)
            self._set_row_visible("half_thickness", True)
            self._on_impedance_changed(
                self._widgets["impedance"].currentText())

    def _on_impedance_changed(self, imp):
        has_wp = bool(self._widgets["workpiece"].text().strip())
        if not has_wp:
            return
        is_esim = (imp == "esim")
        self._set_row_visible("mu_r", not is_esim)
        self._set_row_visible("bh_file", is_esim)
        self._set_row_visible("esim_geometry", is_esim)

    def _on_fem_material_changed(self, idx):
        method = self.val("method")
        if method != "FEM":
            return
        self._set_row_visible("fem_bh_file", idx == 1)

    def _on_validation_changed(self, _text=None):
        method = self.val("method")
        if method == "BEM":
            for key in ("source", "sink"):
                w = self._widgets[key]
                w.setStyleSheet(self._RED if not w.text().strip()
                                else self._NORMAL)
        if callable(self.validationChanged):
            self.validationChanged()

    def is_runnable(self):
        method = self.val("method")
        if method == "BEM":
            return (bool(self._widgets["source"].text().strip())
                    and bool(self._widgets["sink"].text().strip()))
        return True  # FEM has no hard requirements beyond .vol

    def build_command(self, vol_path):
        if not vol_path:
            raise ValueError("No .vol file specified.")
        method = self.val("method")

        if method == "BEM":
            return self._build_bem_command(vol_path)
        else:
            return self._build_fem_command(vol_path)

    def _build_bem_command(self, vol_path):
        src = self.val("source")
        sink = self.val("sink")
        if not src or not sink:
            raise ValueError("Source and Sink blocks are required.")
        cmd = [_PYTHON, calc_script("calc_inductance.py"),
               "--vol", vol_path,
               "--source", src, "--sink", sink,
               "--frequency", self.val("freq"),
               "--coil-sigma", self.val("coil_sigma")]
        fes = self.val("fes_order")
        if fes and fes != "0":
            cmd += ["--fes-order", fes]
        wp = self.val("workpiece")
        if wp:
            imp = self.val("impedance")
            cmd += ["--workpiece", wp,
                    "--impedance-model", imp,
                    "--wp-sigma", self.val("wp_sigma"),
                    "--half-thickness", self.val("half_thickness")]
            if imp == "esim":
                cmd += ["--esim-geometry", self.val("esim_geometry")]
                bh = self.val("bh_file")
                if bh:
                    cmd += ["--bh-file", bh]
            else:
                cmd += ["--mu-r", self.val("mu_r")]
        air = self.val("air")
        if air:
            cmd += ["--air", air,
                    "--msh-output", msh_output(vol_path, "_bem")]
        return cmd

    def _build_fem_command(self, vol_path):
        mat_idx = self._widgets["fem_material"].currentIndex()
        if mat_idx == 0:
            mat_args = ["--material", "custom",
                        "--mu-r", self.val("mu_r")]
        else:
            mat_args = ["--material", "steel"]
            bh = self.val("fem_bh_file")
            if bh:
                mat_args += ["--bh-file", bh]
        cmd = [_PYTHON, calc_script("calc_fem_kelvin.py"),
               "--vol", vol_path,
               "--fes-order", self.val("fes_order"),
               "--frequency", self.val("freq"),
               "--sigma", self.val("wp_sigma"),
               "--impedance", self.val("impedance"),
               "--current", self.val("current"),
               "--a-coil", self.val("a_coil"),
               "--r-wp", self.val("r_wp"),
               "--solver", self.val("solver"),
               "--max-iter", self.val("max_iter"),
               "--msh-output", msh_output(vol_path, "_fem")]
        cmd += mat_args
        return cmd


class IHWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__("Radia - Induction Heating", vol_path,
                         settings_key="ih")
        panel = IHPanel()
        self._set_panel(panel)
        self._restore_settings()
        self._update_run_state()


def main():
    run_app(IHWindow)


if __name__ == "__main__":
    main()
