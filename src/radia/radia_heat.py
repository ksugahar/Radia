"""Radia thermal analysis window (Phase B of the IH pipeline).

Two-stage workflow with the IH window:

  1. Run an IH method (PEEC+BEM / PEEC+FEM Kelvin) in radia_ih.py
     against the EM .vol.  This writes <stem>_qsurf.sol, the surface
     heat flux distribution.
  2. Run radia_heat.py against a SEPARATE workpiece-volume .vol
     (the IH SIBC mesh treats the workpiece as a hole; the thermal
     mesh has it as a real solid).  Provide the qsurf .sol from
     step 1 to pick up the spatial distribution.

The panel exposes a "Heat source" selector with two values:

  - "Uniform q_surf [W/m^2]" : single scalar applied across the whole
    heating face.  Cheap.  Useful for first-cut feasibility runs and
    smoke testing the thermal pipeline without paying for an EM solve.
  - "Spatial q_surf .sol (from IH)" : load the .sol that
    calc_fem_kelvin.py emits.  Preserves the spatial distribution
    (hotspot under the coil, cold zones away from it).

Time integration: backward Euler (default) or Crank-Nicolson.

Standalone launchable -- this module wires up an
:class:`AnalysisWindow` exactly like radia_ih.py, so the
"radia-heat" entry-point script can spawn it without Cubit.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import (
    ModePanel, AnalysisWindow, calc_script, msh_output, json_output,
    run_app, _PYTHON,
)
from PySide6.QtWidgets import QFileDialog


# ============================================================
# Constants
# ============================================================

TITLE = "Thermal"
# We do not require any specific .vol label set here -- the user
# picks --surface-label in the panel.  Just hint the common ones.
REQUIRED_LABELS = []
OPTIONAL_LABELS = ["work_main", "work_skin", "outer", "sibc"]
OPTIONAL_FILES = {}


HEAT_SRC_UNIFORM = "Uniform q_surf [W/m^2]"
HEAT_SRC_SPATIAL = "Spatial q_surf .sol (from IH)"

MESH_TYPE_3D     = "3D volume"
MESH_TYPE_AXISYM = "2D axisymmetric (r, z)"


# Material thermal preset labels (UI-side; calc_heat.py knows the
# numeric values).  Keep the order consistent with calc_heat.py's
# THERMAL_PRESETS so the index round-trip via QSettings is stable.
THERMAL_PRESET_NAMES = [
    "Steel (Fe)", "Aluminum", "Copper", "Stainless 304", "Brass", "Custom",
]
THERMAL_PRESET_TO_CLI = {
    "Steel (Fe)":     "steel",
    "Aluminum":       "aluminum",
    "Copper":         "copper",
    "Stainless 304":  "stainless",
    "Brass":          "brass",
    "Custom":         "custom",
}
THERMAL_PRESET_VALUES = {
    "Steel (Fe)":     {"rho": 7800, "cp": 467, "k": 46.6},
    "Aluminum":       {"rho": 2700, "cp": 900, "k": 237.0},
    "Copper":         {"rho": 8960, "cp": 385, "k": 401.0},
    "Stainless 304":  {"rho": 8000, "cp": 500, "k": 16.0},
    "Brass":          {"rho": 8530, "cp": 380, "k": 109.0},
    "Custom":         {},
}


# ============================================================
# Heat panel
# ============================================================

class HeatPanel(ModePanel):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._vol_mats = None
        self._vol_bnds = None
        self._build_ui()

    # ------- UI construction -------

    def _build_ui(self):
        # Mesh type selection (3D volume vs 2D axisymmetric).
        # This is the routing key: 3D goes to calc_heat.py, axisym
        # goes to calc_heat_axisym.py.  Both consume the same
        # qsurf .sol from the IH solve.
        self._add_section("Mesh type")
        mesh_t = self.add_combo(
            "mesh_type", "Mesh:",
            [MESH_TYPE_3D, MESH_TYPE_AXISYM], default=0)
        mesh_t.currentTextChanged.connect(self._on_mesh_type_changed)
        mesh_t.setToolTip(
            "<b>3D volume</b>: arbitrary 3D workpiece -- runs "
            "calc_heat.py with the volumetric heat equation.<br>"
            "<b>2D axisymmetric</b>: rotationally symmetric "
            "workpiece (cylinder, stepped shaft) meshed in the "
            "(r, z) plane -- runs calc_heat_axisym.py with "
            "(2*pi*r) weighting.  10-100x faster than the "
            "equivalent 3D mesh for typical IH cylinder cases.<br>"
            "Cross-mesh q_surf transfer is phi-averaged in axisym "
            "mode so a slightly non-axisymmetric coil (gapped "
            "torus) still produces a physically sensible q.")
        self.add_spin("n_phi_samples", "phi samples (axisym):", 8, 1, 64)

        # Heat source selection (the solver-switch variable).
        self._add_section("Heat source")
        src = self.add_combo(
            "heat_source", "Source:",
            [HEAT_SRC_UNIFORM, HEAT_SRC_SPATIAL])
        src.currentTextChanged.connect(self._on_heat_source_changed)
        src.setToolTip(
            "<b>Uniform</b>: a single scalar [W/m^2] applied across "
            "the whole heating face.  Cheap, useful for testing.<br>"
            "<b>Spatial .sol</b>: load q_surf saved by an IH solve "
            "(calc_fem_kelvin.py).  Preserves hotspot distribution.")

        self.add_line("q_uniform", "q_surf uniform [W/m^2]:", "1.0e6",
                      placeholder="e.g. 1e6 W/m^2")

        # Spatial-mode inputs.  Visibility toggled by
        # _on_heat_source_changed.
        self._add_section("Spatial q_surf source", key="_sec_spatial")
        self.add_browse(
            "qsurf_sol", "qsurf .sol:",
            filter_str="NGSolve sol (*.sol);;All (*)")
        self.add_browse(
            "em_vol", "EM .vol (auto):",
            filter_str="Netgen volume (*.vol);;All (*)")
        self.add_spin("qsurf_order", "qsurf H1 order:", 1, 1, 5)

        # Workpiece thermal mesh.
        self._add_section("Workpiece thermal mesh")
        wp_w = self.add_browse(
            "wp_vol", "wp .vol:",
            filter_str="Netgen volume (*.vol);;All (*)")
        wp_w.textChanged.connect(self._emit_validation)
        # Surface label entry: a .vol-derived combo would be nicer
        # but parity with the calc CLI argument is the priority.
        self.add_line("surface_label", "Heating surface label:", "outer")

        # Material.
        self._add_section("Material (workpiece thermal)")
        mat = self.add_combo(
            "thermal_material", "Preset:",
            THERMAL_PRESET_NAMES, default=0)  # Steel
        mat.currentTextChanged.connect(self._on_thermal_material_changed)
        mat.setToolTip(
            "Workpiece thermal properties.  Custom: enter rho, cp, k "
            "manually.  Presets are room-temperature values; for "
            "high-temperature accuracy use Custom and dial in the "
            "values from a property table.")
        self.add_line("rho", "rho [kg/m^3]:", "7800")
        self.add_line("cp",  "cp [J/(kg.K)]:", "467")
        self.add_line("k",   "k [W/(m.K)]:", "46.6")

        # Boundary conditions.
        self._add_section("Boundary conditions")
        self.add_line("h_conv", "h_conv [W/(m^2.K)]:", "10")
        self.add_line("t_ext",  "T_ext [degC]:",       "20")
        self.add_line("t_init", "T_initial [degC]:",   "20")

        # Time integration (the second solver-switch).
        self._add_section("Time integration")
        scheme = self.add_combo(
            "time_scheme", "Scheme:",
            ["Backward Euler", "Crank-Nicolson"], default=0)
        scheme.setToolTip(
            "<b>Backward Euler</b> (theta=1, default): unconditionally "
            "stable, first-order accurate in time.  Matches Kubota's "
            "notebook.<br>"
            "<b>Crank-Nicolson</b> (theta=0.5): unconditionally stable, "
            "second-order accurate.  Slightly oscillatory if dt is "
            "large relative to the thermal time-constant.")
        self.add_line("dt",    "dt [s]:",    "0.5")
        self.add_line("t_end", "t_end [s]:", "5.0")

        # Linear solver.
        self._add_section("Linear solver")
        self.add_combo(
            "linear_solver", "Solver:",
            ["sparsecholesky", "umfpack", "pardiso"], default=0)

        # Observation.
        self._add_section("Observation")
        probe = self.add_line("probe_point",
                              "Probe point x,y,z [m]:", "")
        probe.setToolTip(
            "Optional probe point for the T(t) history.  Format: "
            "comma-separated x,y,z in metres.  Example: '0.0305,0,0'.")
        self.add_browse(
            "csv_output", "CSV out:", default="",
            filter_str="CSV (*.csv);;All (*)")

        # Advanced.
        self._add_section("Advanced")
        self.add_spin("fes_order", "H1 order:", 1, 1, 3)
        self.add_browse(
            "vtu_prefix", "VTU prefix:", default="",
            filter_str="(no extension);;All (*)")

        # Initial visibility / preset.
        self._on_mesh_type_changed(mesh_t.currentText())
        self._on_heat_source_changed(src.currentText())
        self._on_thermal_material_changed(mat.currentText())

    def _on_mesh_type_changed(self, name):
        is_axisym = (name == MESH_TYPE_AXISYM)
        # n_phi_samples is meaningful only in axisym mode.
        self._set_row_visible("n_phi_samples", is_axisym)

    # ------- Handlers -------

    def _emit_validation(self, *_):
        cb = getattr(self, "validationChanged", None)
        if callable(cb):
            cb()

    def _on_heat_source_changed(self, name):
        is_uniform = (name == HEAT_SRC_UNIFORM)
        # Toggle spatial section + uniform line.  ``_sec_spatial``
        # is the section header key; collapsing it removes the
        # entire group when the uniform-source mode is active.
        self._set_row_visible("q_uniform", is_uniform)
        self._set_row_visible("_sec_spatial", not is_uniform)
        for key in ("qsurf_sol", "em_vol", "qsurf_order"):
            self._set_row_visible(key, not is_uniform)

    def _on_thermal_material_changed(self, name):
        cli = THERMAL_PRESET_TO_CLI.get(name, "custom")
        is_custom = (cli == "custom")
        preset = THERMAL_PRESET_VALUES.get(name, {})
        for key in ("rho", "cp", "k"):
            w = self._widgets.get(key)
            if w is None:
                continue
            if not is_custom and key in preset:
                w.setText(f"{preset[key]:g}")
            w.setEnabled(is_custom)

    # ------- Validation / command building -------

    def wp_vol_path(self):
        """Used by AnalysisWindow._on_run() to set subprocess work_dir."""
        return self.val("wp_vol") if "wp_vol" in self._widgets else ""

    def is_runnable(self):
        wp = self.val("wp_vol") if "wp_vol" in self._widgets else ""
        if not (wp and os.path.isfile(wp)):
            return False
        if self.val("heat_source") == HEAT_SRC_UNIFORM:
            try:
                float(self.val("q_uniform"))
                return True
            except ValueError:
                return False
        # Spatial mode: need a .sol; em-vol auto-detects.
        sol = self.val("qsurf_sol")
        return bool(sol and os.path.isfile(sol))

    def build_command(self, vol_path):
        # vol_path is the Window-level .vol from the launcher; we
        # ignore it because Heat takes its own --wp-vol.
        wp = self.val("wp_vol")
        if not wp:
            raise ValueError("Workpiece .vol is required for Heat.")
        if not os.path.isfile(wp):
            raise ValueError(f"wp .vol not found: {wp}")

        material_cli = THERMAL_PRESET_TO_CLI.get(
            self.val("thermal_material"), "custom")
        scheme_map = {
            "Backward Euler": "backward-euler",
            "Crank-Nicolson": "crank-nicolson",
        }
        scheme_cli = scheme_map.get(
            self.val("time_scheme"), "backward-euler")

        # Route to the correct calc script based on mesh type.
        is_axisym = (self.val("mesh_type") == MESH_TYPE_AXISYM)
        calc = ("calc_heat_axisym.py" if is_axisym
                else "calc_heat.py")
        cmd = [_PYTHON, calc_script(calc),
               "--wp-vol", wp,
               "--surface-label", self.val("surface_label"),
               "--material", material_cli,
               "--rho", self.val("rho"),
               "--cp",  self.val("cp"),
               "--k",   self.val("k"),
               "--h-conv",   self.val("h_conv"),
               "--t-ext",    self.val("t_ext"),
               "--t-initial",self.val("t_init"),
               "--dt",       self.val("dt"),
               "--t-end",    self.val("t_end"),
               "--time-scheme", scheme_cli,
               "--linear-solver", self.val("linear_solver"),
               "--fes-order", str(self.val("fes_order")),
               "--msh-output", msh_output(wp, "_heat"),
               "--output", json_output(wp, "_heat")]

        # Heat source.
        if self.val("heat_source") == HEAT_SRC_UNIFORM:
            cmd += ["--q-uniform", self.val("q_uniform")]
        else:
            cmd += ["--qsurf-sol", self.val("qsurf_sol"),
                    "--qsurf-order", str(self.val("qsurf_order"))]
            em_vol = self.val("em_vol")
            if em_vol:
                cmd += ["--em-vol", em_vol]
            if is_axisym:
                cmd += ["--n-phi-samples", str(self.val("n_phi_samples"))]

        # Probe + CSV.
        probe = self.val("probe_point").strip()
        if probe:
            cmd += ["--probe-point", probe]
        csv = self.val("csv_output")
        if csv:
            cmd += ["--csv-output", csv]
        vtu = self.val("vtu_prefix")
        if vtu:
            cmd += ["--vtu-prefix", vtu]

        return cmd


# ============================================================
# Heat window
# ============================================================

class HeatWindow(AnalysisWindow):
    def __init__(self, vol_path="", *, prefill=None):
        super().__init__("Radia - Thermal", vol_path,
                         settings_key="heat")
        panel = HeatPanel()
        self._set_panel(panel)
        self._restore_settings()
        # If the launcher passed a .vol, treat it as the wp .vol
        # (the most common chain-from-IH case).  Display the path
        # relative to the working folder if it lives inside it.
        if vol_path and "wp_vol" in panel._widgets and \
                not panel.val("wp_vol"):
            panel._widgets["wp_vol"].setText(self.display_path(vol_path))
        # Apply chain-launch pre-fill (qsurf_sol, em_vol, ...).  The
        # IH window passes these via the CLI when the user clicks
        # "Run thermal..." after a successful IH solve so the user
        # does not have to type the paths in twice.
        if prefill:
            self._apply_prefill(panel, prefill)
        self._update_run_state()

    @staticmethod
    def _apply_prefill(panel, prefill):
        """Set widget values from a dict of {key: value}.

        Auto-flips ``heat_source`` to spatial mode when a qsurf-sol
        is present so the user lands directly on a runnable panel.
        """
        for key, value in prefill.items():
            if not value:
                continue
            w = panel._widgets.get(key)
            if w is None:
                continue
            if hasattr(w, "setText"):
                w.setText(str(value))
            elif hasattr(w, "setCurrentText"):
                w.setCurrentText(str(value))
            elif hasattr(w, "setValue"):
                try:
                    w.setValue(int(value))
                except (TypeError, ValueError):
                    pass
        if prefill.get("qsurf_sol"):
            src = panel._widgets.get("heat_source")
            if src is not None and hasattr(src, "setCurrentText"):
                src.setCurrentText(HEAT_SRC_SPATIAL)


def main():
    """Standalone entry point.

    Usage::

        python radia_heat.py [WP_VOL] [--qsurf-sol PATH]
                              [--em-vol PATH] [--qsurf-order N]

    The optional flags pre-fill the panel; they are how the IH
    window's "Run thermal..." chain button passes its EM-side
    output through.
    """
    import argparse
    parser = argparse.ArgumentParser(
        description="Radia thermal analysis panel.",
        # Don't intercept --help inside the GUI subprocess; let it
        # fall through to argparse's normal help-and-exit so the
        # caller can discover the chain-fill flags.
    )
    parser.add_argument("wp_vol", nargs="?", default="",
                        help="Optional workpiece .vol path (pre-fills "
                             "the wp .vol field).")
    parser.add_argument("--qsurf-sol", default="",
                        help="Pre-fill the qsurf .sol field "
                             "(typically <em-msh stem>_qsurf.sol "
                             "from a calc_fem_kelvin.py run).")
    parser.add_argument("--em-vol", default="",
                        help="Pre-fill the EM .vol field "
                             "(typically <em-msh stem>_fem.vol).")
    parser.add_argument("--qsurf-order", type=int, default=None,
                        help="Pre-fill the qsurf H1 order spinner.")
    args = parser.parse_args()

    prefill = {}
    if args.qsurf_sol:
        prefill["qsurf_sol"] = args.qsurf_sol
    if args.em_vol:
        prefill["em_vol"] = args.em_vol
    if args.qsurf_order is not None:
        prefill["qsurf_order"] = args.qsurf_order

    from PySide6.QtWidgets import QApplication
    from radia_gui_base import apply_panel_base_font
    app = QApplication([sys.argv[0]])
    apply_panel_base_font(app)
    window = HeatWindow(args.wp_vol, prefill=prefill or None)
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
