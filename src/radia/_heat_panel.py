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
# Helpers
# ============================================================

def _parse_vol_bcnames(vol_path):
    """Return the list of boundary-condition names in a Netgen .vol.

    Parses the ``bcnames`` section from a Netgen .vol text file.  Used
    by the heat panel's surface_label auto-detect (avoids importing
    ngsolve just for metadata, saves ~1-2 s of import latency on every
    Browse-keystroke).

    Format (Netgen .vol):

        bcnames
        <count>
        1 first_label
        2 second_label
        ...

    Returns an empty list on any parse failure or when the file has no
    ``bcnames`` section.
    """
    names = []
    with open(vol_path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() == "bcnames":
            try:
                count = int(lines[i + 1].strip())
            except (ValueError, IndexError):
                return []
            for j in range(count):
                parts = lines[i + 2 + j].strip().split(None, 1)
                if len(parts) == 2:
                    names.append(parts[1])
            return names
        i += 1
    return names


def _parse_vol_materials(vol_path):
    """Return the list of volume-material names in a Netgen .vol.

    Parses the ``materials`` section (same line layout as ``bcnames``)
    from a Netgen .vol text file.  Used by the heat panel to warn early
    when the user Browse-selects a multi-material (coil+wp) mesh: the
    thermal step targets the WORKPIECE SOLID only, so a mesh with more
    than one volume material is rejected by calc_heat.py at Run time.
    Reading the text directly avoids the ~1-2 s ngsolve import on every
    Browse-keystroke.

    Format (Netgen .vol)::

        materials
        <count>
        1 first_material
        2 second_material
        ...

    Returns an empty list on any parse failure or when the file has no
    ``materials`` section (older exports) -- the panel then stays silent
    and lets the calc_heat.py-side guard surface the error at Run.
    """
    names = []
    try:
        with open(vol_path, encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except OSError:
        return names
    i = 0
    while i < len(lines):
        if lines[i].strip() == "materials":
            try:
                count = int(lines[i + 1].strip())
            except (ValueError, IndexError):
                return []
            for j in range(count):
                parts = lines[i + 2 + j].strip().split(None, 1)
                if len(parts) == 2:
                    names.append(parts[1])
            return names
        i += 1
    return names


def _derive_em_vol_from_qsurf(qsurf_path):
    """Derive the companion EM .vol from a qsurf .sol path.

    calc_fem_kelvin.py writes ``<stem>_qsurf.sol`` alongside
    ``<stem>_fem.vol``.  Given the .sol, this returns the sibling .vol
    ONLY when the convention matches AND the file exists -- otherwise ""
    (never a silent wrong guess).  Used to auto-fill the em_vol field so
    the user specifies just the qsurf .sol + the workpiece .vol.
    """
    if not qsurf_path:
        return ""
    suffix = "_qsurf.sol"
    if qsurf_path.lower().endswith(suffix):
        cand = qsurf_path[:-len(suffix)] + "_fem.vol"
        if os.path.isfile(cand):
            return cand
    return ""


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
        #
        # The "Mesh type" section header takes a key so callers that
        # hide the mesh_type combo (e.g. IHPanel, where the parent
        # Method dropdown encodes the choice) can also hide the
        # section header — otherwise it becomes an orphan that the
        # panel_qa check flags.
        self._add_section("Mesh type", key="_sec_mesh_type")
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
        qs_w = self.add_browse(
            "qsurf_sol", "qsurf .sol:",
            filter_str="NGSolve sol (*.sol);;All (*)")
        qs_w.setToolTip(
            "Surface heat-density q_surf [W/m^2] from the IH EM solve "
            "(calc_fem_kelvin.py).  Picking it AUTO-FILLS the EM .vol "
            "below from the companion ``<stem>_fem.vol`` -- so you "
            "normally specify just the qsurf .sol + the workpiece .vol.")
        qs_w.textChanged.connect(self._on_qsurf_sol_changed)
        em_w = self.add_browse(
            "em_vol", "EM .vol:",
            filter_str="Netgen volume (*.vol);;All (*)")
        em_w.setToolTip(
            "EM .vol the qsurf .sol is defined on (the .sol is a "
            "coefficient vector with no embedded mesh, so this is "
            "REQUIRED).  AUTO-FILLED from the qsurf .sol's companion "
            "``<stem>_fem.vol`` when you pick the .sol -- override here "
            "only if your EM mesh is named differently.  The calc script "
            "still receives --em-vol explicitly (no silent fallback).")
        em_w.textChanged.connect(self._emit_validation)
        self.add_spin("qsurf_order", "qsurf H1 order:", 1, 1, 5)
        qpa = self.add_check(
            "q_phi_average",
            "Azimuthal-average q_surf (uniform / axisymmetric):", False)
        qpa.setToolTip(
            "<b>3D solver only.</b>  Circumferentially (phi) average the "
            "spatial q_surf into an AXISYMMETRIC heat input on the 3D "
            "workpiece mesh, then solve WITHOUT rotation time-stepping. "
            "This is the 'uniform' (fast-spinning) limit -- the "
            "temperature comes out axisymmetric.<br>"
            "<b>Leave UNCHECKED with Rotation = 0 rpm</b> for the "
            "<b>no-rotation</b> mode: the spatial q_surf is applied as-is, "
            "giving a non-axisymmetric temperature (hot where the coil is "
            "close).<br>"
            "(The 2D axisym mesh mode already phi-averages, so this "
            "checkbox is hidden there.)")

        # Workpiece thermal mesh.
        self._add_section("Workpiece thermal mesh")
        wp_w = self.add_browse(
            "wp_vol", "wp .vol:",
            filter_str="Netgen volume (*.vol);;All (*)")
        wp_w.textChanged.connect(self._emit_validation)
        wp_w.textChanged.connect(self._on_wp_vol_changed_for_bnd)
        # Surface label: editable combo populated from the wp .vol's
        # bcnames when Browse-selected (P1, v4.74.0).  Empty entry
        # means "apply to ALL BND" (P2, calc_heat.py treats empty as
        # the ".*" regex match).  This unifies two friction points:
        # the old default "outer" mismatched keiko/kubota's "sibc"
        # workpieces (they had to retype), and a single-BND workpiece
        # shouldn't need any label at all.
        sl = self.add_combo(
            "surface_label", "Heating surface label:", [""], default=0)
        sl.setEditable(True)
        sl.setToolTip(
            "Boundary label where qsurf and Newton convection are "
            "applied.<br><br>"
            "<b>Leave empty</b> to apply to ALL BND (the common case "
            "for a single-workpiece .vol -- the panel auto-fills the "
            "sole BND label when the .vol has exactly one).<br><br>"
            "Pick a specific label only when the workpiece has "
            "multiple BND sidesets and heating + convection should be "
            "restricted to a subset.")

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
            "Override rho/cp/k", default=False)
        ov.toggled.connect(self._on_override_toggled)
        ov.setToolTip(
            "Use preset values as a starting point and edit rho/cp/k "
            "below.  When OFF, rho/cp/k are locked to the preset's "
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
        #
        # v4.78.0: QDoubleSpinBox replaces QLineEdit so the user can't
        # type non-numeric or accidentally negative values; suffix " rpm"
        # makes the unit unmistakable.  Sign convention: positive =
        # CCW viewed from the +axis end (right-hand rule), enforced by
        # the spinbox lower bound 0.0 (a negative rpm is physically the
        # same as positive rpm around the OPPOSITE axis, so we narrow
        # the input space to one canonical form).
        from radia_gui_base import _NoWheelDoubleSpinBox as QDoubleSpinBox
        rrow = QDoubleSpinBox()
        rrow.setRange(0.0, 1e6)
        rrow.setDecimals(2)
        rrow.setSingleStep(10.0)
        rrow.setSuffix(" rpm")
        rrow.setValue(0.0)
        self._form.addRow("Rotation (0 = stationary):", rrow)
        self._widgets["rotation_rpm"] = rrow
        self._row_indices["rotation_rpm"] = self._form.rowCount() - 1
        rrow.setToolTip(
            "<b>3D solver</b>: positive rpm makes the workpiece body "
            "spin around the chosen axis -- q_surf is re-sampled on "
            "the rotated body each timestep.<br>"
            "<b>2D axisym solver</b>: rotation is implicit (the "
            "workpiece is rotation-symmetric by construction); the "
            "value is recorded as metadata only -- the answer does "
            "NOT depend on rpm in axisym mode.<br>"
            "<b>Uniform source</b>: rotation has no effect (q_surf is "
            "constant); the field stays editable but is ignored.")

        # Rotation axis -- v4.78.0 (was hardcoded to z pre-v4.78.0).
        # Horizontal-axis workpieces (billet along x) silently got the
        # wrong physics under the old z-only assumption.
        axis_combo = self.add_combo(
            "rotation_axis", "Rotation axis:",
            ["z", "x", "y"], default=0)
        axis_combo.setToolTip(
            "Workpiece rotation axis (default z).  Positive rpm = CCW "
            "viewed from the positive end of this axis (right-hand "
            "rule).<br><br>"
            "Verify the axis matches your wp_vol's geometry: print "
            "the bounding-box extents in the BND log and pick the axis "
            "the workpiece is symmetric around.<br><br>"
            "Pre-v4.78.0 the axis was hardcoded to z; horizontal-axis "
            "workpieces silently rotated around the wrong axis.")

        # Boundary conditions.
        self._add_section("Boundary conditions")
        self.add_line("h_conv", "h_conv [W/(m^2.K)]:", "10")
        self.add_line("t_ext",  "T_ext [degC]:",       "20")
        self.add_line("emissivity", "Emissivity [0..1]:", "0")
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
        self._update_q_phi_avg_visible()

    def _update_q_phi_avg_visible(self):
        """The 3D uniform / phi-average checkbox is meaningful only for a
        3D mesh + spatial source: the axisym mesh already phi-averages,
        and a uniform-constant q is already azimuthally uniform."""
        is_axisym = (self.val("mesh_type") == MESH_TYPE_AXISYM)
        is_uniform = (self.val("heat_source") == HEAT_SRC_UNIFORM)
        self._set_row_visible(
            "q_phi_average", (not is_axisym) and (not is_uniform))

    # ------- Handlers -------

    def _emit_validation(self, *_):
        cb = getattr(self, "validationChanged", None)
        if callable(cb):
            cb()

    def _on_wp_vol_changed_for_bnd(self, _path):
        """Auto-populate the surface_label combo from the wp .vol BNDs.

        P1 of the v4.74.0 thermal-panel UX improvement: when the user
        Browse-selects a workpiece .vol, read its bcnames and either
        (a) auto-fill ``surface_label`` if there is exactly one BND
        label, or (b) populate the combo dropdown with the available
        labels (still editable so the user can type a regex).

        The auto-fill case fixes the historical friction where the
        panel default ``'outer'`` didn't match keiko/kubota's actual
        BND name (``'sibc'``) -- they had to manually retype to match
        the .vol.

        On any failure (file missing, parse error, ngsolve not yet
        loaded), this is silent: the combo keeps its current state
        and the calc_heat.py-side validation will surface a clear
        error at Run time.  Empty entry (the first item) always
        remains so the user can intentionally select "all BND" -- P2
        in calc_heat.py interprets that as ``.*``.
        """
        # Resolve the absolute path via the AnalysisWindow's helper
        # (same one ``val()`` uses) so display-relative paths work.
        try:
            vol_path = self.val("wp_vol")
        except Exception:
            return
        if not vol_path or not os.path.isfile(vol_path):
            return

        # Workpiece-only early warning: the Thermal step targets the
        # WORKPIECE SOLID only.  If the .vol has >1 volume material
        # (e.g. a coil+workpiece EM mesh), calc_heat.py rejects it at
        # Run -- flag it now on the wp_vol field so the user sees it on
        # Browse, not after a Run round-trip.  (wp_vol has no base
        # tooltip, so clearing to "" on a clean single-material mesh is
        # safe.)
        wp_widget = self._widgets.get("wp_vol")
        if wp_widget is not None:
            mats = sorted(set(_parse_vol_materials(vol_path)))
            if len(mats) > 1:
                wp_widget.setToolTip(
                    "WARNING: this .vol has {} volume materials {}.  The "
                    "Thermal step targets the WORKPIECE SOLID only and "
                    "will reject a multi-material (coil+workpiece) mesh. "
                    "Export a workpiece-only volume mesh (a single solid) "
                    "for Thermal.".format(len(mats), mats))
            else:
                wp_widget.setToolTip("")

        try:
            # Parse .vol text directly so we don't pay the ngsolve
            # import cost (~1-2s) on every keystroke -- the bcnames
            # section is a small, well-defined block near the top
            # of the file.
            bnds = _parse_vol_bcnames(vol_path)
        except Exception:
            return
        if not bnds:
            return
        # Capture current selection so we can preserve user intent
        # (e.g. they already typed "top" before pasting the .vol path).
        combo = self._widgets["surface_label"]
        prior = combo.currentText().strip()
        # Repopulate: empty first ("apply to ALL BND"), then the
        # sorted unique BND names.
        items = [""] + sorted(set(bnds))
        combo.blockSignals(True)
        try:
            combo.clear()
            combo.addItems(items)
            if prior and prior in items:
                # User typed a valid label before -- keep it.
                combo.setCurrentIndex(items.index(prior))
            elif len(items) == 2:
                # Exactly one BND label (items = ["", name]); auto-
                # fill it.  This is the keiko/kubota single-workpiece
                # case where naming the sole BND is friction.
                combo.setCurrentIndex(1)
            else:
                # Multiple BNDs -- leave empty (apply to ALL by P2).
                combo.setCurrentIndex(0)
        finally:
            combo.blockSignals(False)
        # Re-emit validation in case the new selection changes the
        # is_runnable state (unlikely but free).
        self._emit_validation()

    def _on_qsurf_sol_changed(self, _text):
        """Auto-fill em_vol from the qsurf .sol's companion
        ``<stem>_fem.vol`` so the user specifies just qsurf .sol +
        workpiece .vol.

        Never clobbers a user-supplied em_vol: it only writes when
        em_vol is empty or still holds the last auto value (tracked in
        ``self._em_vol_auto``).  When the new .sol has no derivable
        companion, a stale auto value is CLEARED so the user browses the
        correct EM .vol (fail-fast, no stale guess).  The calc script
        still receives --em-vol explicitly (calc-side No-Fallback intact).
        """
        em_w = self._widgets.get("em_vol")
        if em_w is None:
            return
        cand = _derive_em_vol_from_qsurf(self.val("qsurf_sol"))
        cur = self.val("em_vol")
        prev_auto = getattr(self, "_em_vol_auto", "")
        if cur in ("", prev_auto):
            if cand:
                win = getattr(self, "_analysis_window", None)
                em_w.setText(win.display_path(cand) if win is not None
                             else cand)
                self._em_vol_auto = cand
            elif cur == prev_auto and prev_auto:
                # New .sol has no companion AND em_vol still held the old
                # auto value -> clear the stale path so the user picks the
                # correct EM .vol rather than silently reusing the wrong one.
                em_w.setText("")
                self._em_vol_auto = ""
        self._emit_validation()

    def _on_heat_source_changed(self, name):
        is_uniform = (name == HEAT_SRC_UNIFORM)
        # Toggle spatial section + uniform line.  ``_sec_spatial``
        # is the section header key; collapsing it removes the
        # entire group when the uniform-source mode is active.
        self._set_row_visible("q_uniform", is_uniform)
        self._set_row_visible("_sec_spatial", not is_uniform)
        for key in ("qsurf_sol", "em_vol", "qsurf_order"):
            self._set_row_visible(key, not is_uniform)
        self._update_q_phi_avg_visible()
        # v4.78.0: rotation_rpm / rotation_axis are no-ops when source
        # is Uniform (a constant q_surf is rotation-invariant -- the
        # 3D solver's q_resample callback is None on that path).  Grey
        # out the widgets so the user sees that typing a value will be
        # ignored, instead of finding out only by reading the log.
        rpm_w = self._widgets.get("rotation_rpm")
        axis_w = self._widgets.get("rotation_axis")
        if rpm_w is not None:
            rpm_w.setEnabled(not is_uniform)
            if is_uniform:
                rpm_w.setToolTip(
                    "Disabled: Uniform q_surf is rotation-invariant. "
                    "Pick Source = 'Spatial q_surf .sol' to use rotation.")
        if axis_w is not None:
            axis_w.setEnabled(not is_uniform)

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
               "--emissivity", self.val("emissivity"),
               "--t-initial",self.val("t_init"),
               "--dt",       self.val("dt"),
               "--t-end",    self.val("t_end"),
               "--time-scheme", scheme_cli,
               "--linear-solver", self.val("linear_solver"),
               "--fes-order", str(self.val("fes_order")),
               "--rotation-rpm", str(self.val("rotation_rpm")),
               "--rotation-axis", str(self.val("rotation_axis")),
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
                    "Spatial qsurf mode requires a qsurf .sol AND its "
                    "companion EM .vol.  The EM .vol auto-fills from the "
                    ".sol's <stem>_fem.vol sibling when present -- if it "
                    "stayed empty, browse the EM .vol manually (.sol is a "
                    "coefficient vector with no embedded mesh).")
            cmd += ["--qsurf-sol",  sol,
                    "--em-vol",     em_vol,
                    "--qsurf-order", str(self.val("qsurf_order"))]
            if is_axisym:
                cmd += ["--n-phi-samples", str(self.val("n_phi_samples"))]
            elif self._widgets["q_phi_average"].isChecked():
                # 3D uniform / phi-average mode: axisymmetric q on the 3D
                # mesh, no rotation time-stepping.  isChecked() -- NOT
                # self.val() which returns "0"/"1" (both truthy).
                cmd += ["--q-phi-average"]

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
