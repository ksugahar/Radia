"""HeatPanel sub-widget for the radia-ih Thermal method.

Internal implementation detail of ``radia.radia_ih.IHWindow``.  This
panel is embedded as a section that becomes visible when the user
selects ``Method = "Thermal"`` in the IH window.

Promoted from the (now-removed) ``radia.radia_heat`` module in
v4.62.0; the pre-v4.59.0 standalone ``radia-heat`` window has been
retired in favour of the integrated Method-dropdown UX.

This module owns the UI surface of the Thermal method:
  - Mesh type (3D volume vs 2D axisymmetric) -> calc_heat[_axisym].py
  - Heat source (Uniform vs Spatial qsurf .sol)
  - Workpiece thermal mesh / material / convection BC
  - Time integration scheme + dt / t_end
  - Workpiece rotation (v4.58.0+)
  - Probe / CSV / VTU outputs

It does NOT own the actual computation -- that's calc_heat.py /
calc_heat_axisym.py (Layer 4 headless subprocess), invoked via
``build_command()``.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from radia_gui_base import (
    ModePanel, calc_script, msh_output, json_output, _PYTHON,
)


# ============================================================
# Constants
# ============================================================

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
        em_w = self.add_browse(
            "em_vol", "EM .vol:",
            filter_str="Netgen volume (*.vol);;All (*)")
        em_w.setToolTip(
            "EM .vol that the qsurf .sol was saved against.  "
            "REQUIRED -- NGSolve .sol is a coefficient vector only "
            "(no embedded mesh), so the EM .vol must be supplied "
            "explicitly.  Typically the ``<stem>_fem.vol`` that "
            "calc_fem_kelvin.py writes next to ``<stem>_qsurf.sol``.  "
            "Auto-detection was removed 2026-05-20.")
        em_w.textChanged.connect(self._emit_validation)
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
        ov = self.add_check(
            "override_kcprho",
            "Override rho/cp/k (use preset values as starting point, "
            "edit below)", default=False)
        ov.toggled.connect(self._on_override_toggled)
        ov.setToolTip(
            "When OFF, rho/cp/k are locked to the preset's "
            "room-temperature values.  When ON, you can edit any of "
            "the three values directly -- the preset is still emitted "
            "to --material so the JSON output keeps a human-readable "
            "label, and --rho/--cp/--k act as overrides "
            "(calc_heat _resolve_material).")
        self.add_line("rho", "rho [kg/m^3]:", "7800")
        self.add_line("cp",  "cp [J/(kg.K)]:", "467")
        self.add_line("k",   "k [W/(m.K)]:", "46.6")

        # Workpiece rotation.  Since v4.58.0 the 3D solver actually
        # spins the body: q_surf is re-projected each timestep at the
        # workpiece's instantaneous angle around the z axis (mesh /
        # FES / mass / stiffness are held fixed, only the LinearForm
        # RHS reassembles).  For the axisym solver rotation is
        # implicit in the axisymmetric assumption.  Use 0 for a
        # stationary "frozen at one azimuthal configuration" answer
        # (e.g. quick feasibility runs).
        rrow = self.add_line(
            "rotation_rpm",
            "Rotation [rpm] (0 = stationary):", "0")
        rrow.setToolTip(
            "<b>3D solver</b>: positive rpm makes the workpiece body "
            "spin around the z axis -- q_surf is re-sampled on the "
            "rotated body each timestep (uniform mode is unaffected; "
            "only spatial qsurf benefits).<br>"
            "<b>2D axisym solver</b>: rotation is implicit (the "
            "workpiece is rotation-symmetric by construction); the "
            "value is recorded as metadata.")

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
        override = self._is_override_kcprho()
        for key in ("rho", "cp", "k"):
            w = self._widgets.get(key)
            if w is None:
                continue
            # When the preset changes we re-seed the LineEdits with
            # the preset's room-T values, both when locked (showing the
            # values the solver will use) and when Override is on
            # (treating preset as the starting point for the override).
            if not is_custom and key in preset:
                w.setText(f"{preset[key]:g}")
            w.setEnabled(is_custom or override)

    def _is_override_kcprho(self):
        w = self._widgets.get("override_kcprho")
        return bool(w and w.isChecked())

    def _on_override_toggled(self, checked):
        # Re-apply enablement to rho/cp/k without reseeding values: the
        # user wants to keep whatever they've typed when toggling
        # override on/off.
        cli = THERMAL_PRESET_TO_CLI.get(
            self.val("thermal_material"), "custom")
        is_custom = (cli == "custom")
        for key in ("rho", "cp", "k"):
            w = self._widgets.get(key)
            if w is not None:
                w.setEnabled(is_custom or bool(checked))

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
        # Spatial mode: need BOTH a .sol AND its companion .vol -- the
        # .sol is a coefficient vector only, no embedded mesh.  Both
        # paths must point at existing files; otherwise the Run button
        # stays disabled.
        sol = self.val("qsurf_sol")
        emv = self.val("em_vol")
        return bool(sol and os.path.isfile(sol)
                    and emv and os.path.isfile(emv))

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
               "--h-conv",   self.val("h_conv"),
               "--t-ext",    self.val("t_ext"),
               "--t-initial",self.val("t_init"),
               "--dt",       self.val("dt"),
               "--t-end",    self.val("t_end"),
               "--time-scheme", scheme_cli,
               "--linear-solver", self.val("linear_solver"),
               "--fes-order", str(self.val("fes_order")),
               "--rotation-rpm", self.val("rotation_rpm"),
               "--msh-output", msh_output(wp, "_heat"),
               "--output", json_output(wp, "_heat")]

        # Material overrides.  Always pass --rho/--cp/--k when the
        # selected preset is "custom" (then they are required) or when
        # the user has ticked the Override checkbox (then they shadow
        # the preset values via calc_heat._resolve_material).  When
        # neither, we omit them so the JSON output's "rho_kg_m3" etc.
        # report the preset's canonical room-T values rather than a
        # silently echoed copy.
        if material_cli == "custom" or self._is_override_kcprho():
            cmd += ["--rho", self.val("rho"),
                    "--cp",  self.val("cp"),
                    "--k",   self.val("k")]

        # Heat source.
        if self.val("heat_source") == HEAT_SRC_UNIFORM:
            cmd += ["--q-uniform", self.val("q_uniform")]
        else:
            sol = self.val("qsurf_sol")
            em_vol = self.val("em_vol")
            if not (sol and em_vol):
                raise ValueError(
                    "Spatial qsurf mode requires BOTH a qsurf .sol AND "
                    "its companion EM .vol.  .sol files are coefficient "
                    "vectors only -- the .vol carries the mesh.")
            cmd += ["--qsurf-sol",  sol,
                    "--em-vol",     em_vol,
                    "--qsurf-order", str(self.val("qsurf_order"))]
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
