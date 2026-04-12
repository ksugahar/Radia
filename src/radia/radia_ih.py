"""
Radia IH (Induction Heating) analysis window.

Modes:
  BEM          -- EFIE source/sink coil + optional workpiece SIBC
  BEM-SIBC (WP) -- Analytical coil + workpiece surface BEM-SIBC
  FEM          -- Kelvin + SIBC/ESIM

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
REQUIRED_LABELS = ["coil", "source", "sink"]
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
            "method", "Method:", ["BEM", "BEM-SIBC (WP)", "FEM"])
        self._method_combo.currentTextChanged.connect(self._on_method_changed)

        self.add_spin("fes_order", "FES order:", 0, 0, 5)

        # Source / Sink / Coil are NOT GUI inputs. The .jou file is
        # contractually required to use sideset/block names "source",
        # "sink", "coil". calc_inductance.py picks them up by default.
        # Students verify the .jou — after .vol export, no label input
        # is needed at the GUI or calc-script level.

        # BEM saddle-point solver (BEM mode only). Default LU.
        self.add_combo("bem_solver", "BEM solver:",
                       ["LU (direct)", "MINRES", "GMRES"])

        # Coil parameters (always visible)
        self.add_line("freq", "Frequency [Hz]:", "50000")
        self.add_line("coil_sigma", "Coil sigma [S/m]:", "5.8e7")

        # BEM-SIBC (WP) specific: analytical coil geometry
        self.add_line("wp_coil_radius", "Coil radius [m]:", "0.030")
        self.add_line("wp_coil_current", "Coil current [A]:", "1.0")
        self.add_line("wp_gap_deg", "Coil gap [deg]:", "5")
        wp_label_w = self.add_line("wp_bnd_label", "WP boundary label:",
                                   "wp_surface",
                                   placeholder="e.g. wp_surface")
        wp_label_w.textChanged.connect(self._on_validation_changed)

        # Workpiece coupling mode (BEM mode only):
        #   off  = workpiece NOT in DOF (coil-only inductance)
        #   SIBC = include workpiece via linear surface impedance
        #          with per-panel curvature (mesh-driven local R)
        #   ESIM = include workpiece via nonlinear 1D cell problem
        #          with per-panel curvature
        wp_combo = self.add_combo("workpiece_mode", "Workpiece:",
                                  ["off", "SIBC", "ESIM"])
        wp_combo.currentTextChanged.connect(self._on_workpiece_changed)

        # Impedance model combo: kept for FEM and BEM-SIBC (WP) modes
        # which still need an explicit choice. In BEM coil mode the
        # workpiece_mode combo above already determines this, so the
        # impedance row is hidden. Per-panel curvature is now the
        # default (was selectable, simplified 2026-04-12) so the
        # user-facing choice is reduced to "sibc" / "esim".
        imp_combo = self.add_combo("impedance", "Impedance model:",
                                   ["sibc", "esim"])
        imp_combo.currentTextChanged.connect(self._on_impedance_changed)

        self.add_line("wp_sigma", "WP sigma [S/m]:", "2e6")
        self.add_line("half_thickness", "Half thickness [m]:", "0.005")
        self.add_line("mu_r", "mu_r:", "100")
        self.add_browse("bh_file", "BH file:",
                        filter_str="Text files (*.txt *.csv);;All (*)")
        # ESIM 1D cell-problem coordinate system:
        #   "local_curvature": cylindrical Bessel (radial coordinate)
        #   "none": flat slab cosh/sinh (no curvature)
        self.add_combo("esim_geometry", "ESIM geometry:",
                       ["local_curvature", "none"])
        self.add_combo("wp_material", "WP material:",
                       ["steel", "copper", "aluminum"])

        # FEM-specific
        self.add_line("current", "Current [A]:", "1.0")
        self.add_line("a_coil", "Coil radius [m]:", "0.003")
        self.add_line("r_wp", "Workpiece radius [m]:", "0.010")
        self.add_combo("solver", "Solver:",
                       ["pardiso", "bddc", "iccg", "ams"])
        self.add_spin("max_iter", "Max iterations:", 15, 1, 200)

        # BEM-specific: air domain field post-processing on/off.
        # When "on", calc_inductance.py runs B-field post-processing
        # in the air block after the BEM solve.
        self.add_combo("air_mode", "Air field calc:", ["off", "on"])

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
        is_bem_coil = (method == "BEM")
        is_bem_wp = (method == "BEM-SIBC (WP)")
        is_fem = (method == "FEM")

        # BEM coil row: air on/off + bem solver + workpiece mode
        self._set_row_visible("air_mode", is_bem_coil)
        self._set_row_visible("bem_solver", is_bem_coil)

        # BEM coil-specific
        self._set_row_visible("coil_sigma", is_bem_coil)

        # BEM-SIBC (WP) fields (analytical coil, no source/sink)
        for key in ("wp_coil_radius", "wp_coil_current", "wp_gap_deg",
                     "wp_bnd_label"):
            self._set_row_visible(key, is_bem_wp)

        # FEM fields
        for key in ("current", "a_coil", "r_wp", "solver",
                     "max_iter", "fem_material", "fem_bh_file"):
            self._set_row_visible(key, is_fem)

        # Workpiece combo line: BEM coil only (WP mode always has workpiece)
        self._set_row_visible("workpiece_mode", is_bem_coil)

        # Impedance options: only "sibc" and "esim", with per-panel
        # curvature implied. The combo is rebuilt unconditionally so
        # that the previous selection is preserved across method
        # switches when both items are present.
        combo = self._widgets["impedance"]
        prev = combo.currentText()
        combo.clear()
        combo.addItems(["sibc", "esim"])
        idx = combo.findText(prev)
        if idx >= 0:
            combo.setCurrentIndex(idx)

        # BEM-SIBC (WP) always shows workpiece params
        if is_bem_wp:
            for key in ("impedance", "wp_sigma", "half_thickness",
                         "esim_geometry", "wp_material"):
                self._set_row_visible(key, True)
            self._set_row_visible("fes_order", True)
            self._on_impedance_changed(
                self._widgets["impedance"].currentText())
        else:
            self._set_row_visible("wp_material", False)
            self._on_workpiece_changed()

        self._on_validation_changed()

    def _on_workpiece_changed(self, _text=None):
        wp_mode = self.val("workpiece_mode")
        method = self.val("method")
        # workpiece_mode combo is only shown in BEM coil mode; in FEM /
        # BEM-SIBC (WP) modes the workpiece params must stay hidden
        # regardless of what value the combo holds from a previous mode.
        is_bem_coil = (method == "BEM")
        has_wp = is_bem_coil and (wp_mode != "off")

        # Hide all wp params first. The impedance combo is also tied
        # to the workpiece — when wp is off (or in FEM mode where
        # workpiece_mode is hidden) it should disappear so the user
        # can't pick a model that does nothing. In BEM coil mode the
        # impedance row is suppressed entirely because workpiece_mode
        # already determines the impedance choice (SIBC or ESIM).
        for key in ("impedance", "wp_sigma", "half_thickness",
                     "mu_r", "bh_file", "esim_geometry"):
            self._set_row_visible(key, False)

        if has_wp:
            self._set_row_visible("wp_sigma", True)
            self._set_row_visible("half_thickness", True)
            # In BEM coil mode the workpiece_mode combo IS the
            # impedance choice, so the impedance row stays hidden.
            # Map combo to impedance model:
            #   SIBC -> sibc (linear surface impedance, per-panel R)
            #   ESIM -> esim (nonlinear 1D cell problem, per-panel R)
            is_esim = (wp_mode == "ESIM")
            self._set_row_visible("mu_r", not is_esim)
            self._set_row_visible("bh_file", is_esim)
            self._set_row_visible("esim_geometry", is_esim)

    def _on_impedance_changed(self, imp):
        method = self.val("method")
        is_wp_mode = (method == "BEM-SIBC (WP)")
        wp_mode = self.val("workpiece_mode") if "workpiece_mode" in self._widgets else "off"
        has_wp = is_wp_mode or (wp_mode != "off")
        if not has_wp:
            return
        is_esim = (imp == "esim")
        self._set_row_visible("mu_r", not is_esim and not is_wp_mode)
        self._set_row_visible("bh_file", is_esim)
        self._set_row_visible("esim_geometry", is_esim)

    def _on_fem_material_changed(self, idx):
        method = self.val("method")
        if method != "FEM":
            return
        self._set_row_visible("fem_bh_file", idx == 1)

    def _on_validation_changed(self, _text=None):
        cb = getattr(self, 'validationChanged', None)
        if callable(cb):
            cb()

    def is_runnable(self):
        # All methods rely on .jou naming convention; no GUI label inputs.
        return True

    def build_command(self, vol_path):
        if not vol_path:
            raise ValueError("No .vol file specified.")
        method = self.val("method")

        if method == "BEM":
            return self._build_bem_command(vol_path)
        elif method == "BEM-SIBC (WP)":
            return self._build_bem_wp_command(vol_path)
        else:
            return self._build_fem_command(vol_path)

    def _build_bem_command(self, vol_path):
        # Source/sink/coil labels follow the .jou naming convention
        # ("source"/"sink"/"coil") and are picked up by calc_inductance.py
        # defaults. No label arguments are passed from the GUI.
        cmd = [_PYTHON, calc_script("calc_inductance.py"),
               "--vol", vol_path,
               "--frequency", self.val("freq"),
               "--msh-output", msh_output(vol_path, "_bem")]
        coil_sigma = self.val("coil_sigma")
        if coil_sigma:
            cmd += ["--coil-sigma", coil_sigma]
        fes = self.val("fes_order")
        if fes and fes != "0":
            cmd += ["--fes-order", fes]

        # BEM saddle-point solver: combo display name -> --solver value
        bem_solver_map = {
            "LU (direct)": "lu",
            "MINRES": "minres",
            "GMRES": "gmres",
        }
        cmd += ["--solver", bem_solver_map.get(self.val("bem_solver"), "lu")]

        # Workpiece coupling: combo -> --workpiece + impedance
        # The workpiece sideset is named "sibc" by .jou convention
        # (see ih_bem_sample.jou). The label is passed as the
        # boundary name prefix to _extract_panels_from_mesh.
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
            # Per-panel local curvature is now the default for both
            # SIBC and ESIM (was a user-selectable combo, simplified
            # 2026-04-12 — the global half_thickness is only used as
            # a fallback when the extractor cannot recover the local
            # radius, e.g. on degenerate sliver triangles).
            cmd += ["--use-local-curvature"]

        # Air field post-processing on/off
        if self.val("air_mode") == "on":
            cmd += ["--field-air"]

        return cmd

    def _build_bem_wp_command(self, vol_path):
        cmd = [_PYTHON, calc_script("calc_heating_bem.py"),
               "--vol", vol_path,
               "--coil-radius", self.val("wp_coil_radius"),
               "--coil-current", self.val("wp_coil_current"),
               "--gap-deg", self.val("wp_gap_deg"),
               "--frequency", self.val("freq"),
               "--sigma", self.val("wp_sigma"),
               "--material", self.val("wp_material"),
               "--half-thickness", self.val("half_thickness"),
               "--wp-label", self.val("wp_bnd_label"),
               "--msh-output", msh_output(vol_path, "_bem_wp")]
        fes = self.val("fes_order")
        if fes and fes != "0":
            cmd += ["--h1-order", fes]
        imp = self.val("impedance")
        if imp == "esim":
            cmd += ["--esim-geometry", self.val("esim_geometry")]
            bh = self.val("bh_file")
            if bh:
                cmd += ["--bh-file", bh]
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
