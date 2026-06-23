"""Radia IH (Induction Heating) analysis window.

Six methods (Layer-4 dispatch table):

  PEEC inductance (coil only, STEP)
    -> calc_inductance.py --coil-solver peec (no --vol)         (~30 s, L+R_coil)

  BEM-A inductance (coil only, .vol)
    -> calc_inductance.py --coil-solver bem-a (no --vol)        (~25 s, L+R_coil)

  PEEC + BEM weak coupling (workpiece)
    -> calc_inductance.py --coil-solver peec --vol <wp.vol>     (~3 min, L+ΔL+P_wp)

  BEM-A + BEM weak coupling (workpiece)
    -> calc_inductance.py --coil-solver bem-a --vol <wp.vol>    (~3-5 min, L+ΔL+P_wp)

  PEEC coil + FEM wp (SIBC) + Kelvin
    -> calc_fem_kelvin.py --formulation total --peec-step ...   (~4-8 min)

  Full simulation (FEM A-V + wp SIBC + Kelvin)
    -> calc_fem_coilmesh.py                                     (~1-7 min)

Coil topology: gapped torus (real IH has physical port terminations;
closed-torus is unsupported).  Workpiece-coupled modes use weak
coupling: Telegen φ·B back-reaction at the port, coil current
distribution FIXED.  ``calc_peec_bem.py`` was unified into
``calc_inductance.py`` in v4.25.0 (2026-05).
"""
import math
import os
import sys

# Ensure samples directory is reachable from Cubit panels even when
# run from a site-packages install (no editable symlink).
TITLE = "Induction Heating"
REQUIRED_LABELS = []
OPTIONAL_LABELS = ["coil", "air", "kelvin", "kelvin_int", "kelvin_ext",
                    "sibc", "source", "sink", "coil_surface"]
OPTIONAL_FILES = {}

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

try:
    from radia_gui_base import (
        ModePanel, AnalysisWindow, calc_script, msh_output, json_output,
        run_app, _PYTHON,
    )
    from PySide6.QtWidgets import (QLabel, QCheckBox, QGroupBox, QVBoxLayout,
                                     QFormLayout, QWidget)
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor
except ImportError as e:
    # Most likely cause: ``pip install radia[cubit]`` was used without
    # the ``[gui]`` extra, so PySide6 is missing.  The Cubit Solve menu
    # launches this panel via subprocess and would otherwise silently
    # fail with no visible error.  Print a clear install hint to stderr.
    sys.stderr.write(
        "Radia IH panel requires PySide6 but it could not be imported:\n"
        "  {}\n\n"
        "Install with:\n"
        "  pip install --upgrade 'radia[cubit,gui]'\n".format(e))
    sys.exit(1)


# ============================================================
# Constants: method identifiers + labels
# ============================================================

METHOD_PEEC_IND         = "PEEC inductance (coil only, STEP)"
METHOD_BEMA_IND         = "BEM-A inductance (coil only, .vol)"
METHOD_THERMAL_3D_STATIC    = "Thermal: 3D static (no rotation)"
METHOD_THERMAL_3D_ROTATING  = "Thermal: 3D + rotation (q_surf re-sampled per step)"
METHOD_THERMAL_AXISYM       = "Thermal: 2D axisymmetric (rotation implicit)"
THERMAL_METHODS = frozenset({
    METHOD_THERMAL_3D_STATIC,
    METHOD_THERMAL_3D_ROTATING,
    METHOD_THERMAL_AXISYM,
})
METHOD_PEEC_BEM         = "PEEC + BEM weak coupling (workpiece)"
METHOD_BEMA_BEM         = "BEM-A + BEM weak coupling (workpiece)"
METHOD_PEEC_FEM_KELVIN  = "PEEC coil + FEM wp (SIBC) + Kelvin"
METHOD_FEM_FULL         = "Full simulation (FEM A-V + wp SIBC + Kelvin)"

METHOD_TOOLTIP = {
    METHOD_PEEC_IND: (
        "<b>Coil L and R only</b> (~30 s, vacuum).<br>"
        "<ul>"
        "<li>Coil as PEEC filament (no volume mesh, STEP input)</li>"
        "<li>No workpiece, no BEM, no FEM, no .vol needed</li>"
        "<li>Outputs L_coil [nH], R_coil [mOhm] at the given frequency</li>"
        "</ul>"
        "Requires <b>STEP</b> file only.  Use for fast coil impedance "
        "before adding a workpiece."
    ),
    METHOD_BEMA_IND: (
        "<b>Coil L and R only via BEM-A</b> (~25 s, vacuum, surface RWG).<br>"
        "<ul>"
        "<li>Coil as Weggler stabilized EFIE saddle on HDivSurface RWG</li>"
        "<li>n_peri-free (resolves rect-corner current crowding)</li>"
        "<li>Slightly higher L than PEEC (~3% on rect cross-section)</li>"
        "<li>Pre-meshed surface .vol input "
        "(Cubit / Netgen; <b>source</b>/<b>sink</b> sidesets on caps)</li>"
        "</ul>"
        "Requires a coil <b>.vol</b> file only.  Use when PEEC perimeter "
        "filaments cannot resolve the cross-section geometry."
    ),
    METHOD_PEEC_BEM: (
        "<b>Workpiece P + back-reaction Δ L (PEEC coil)</b> "
        "(~3 min, ±5%).<br>"
        "<ul>"
        "<li>Coil as PEEC filament (no volume mesh, STEP input)</li>"
        "<li>Workpiece as BEM-SIBC (surface charge, thin-skin limit)</li>"
        "<li>Weak coupling: Telegen φ·B back-reaction at port "
        "(coil J fixed)</li>"
        "<li>L_total = L_coil_vacuum + Δ L_telegen</li>"
        "</ul>"
        "Requires <b>sibc</b> sideset in .vol + <b>STEP</b> coil."
    ),
    METHOD_BEMA_BEM: (
        "<b>Workpiece P + back-reaction Δ L (BEM-A coil)</b> "
        "(~3-5 min, surface RWG).<br>"
        "<ul>"
        "<li>Coil as BEM-A surface RWG (Weggler EFIE saddle)</li>"
        "<li>Workpiece as BEM-SIBC (surface charge, thin-skin limit)</li>"
        "<li>Weak coupling: Telegen φ·B back-reaction at port</li>"
        "<li>Cross-method check vs PEEC+BEM: agrees ~5% on P_wp</li>"
        "</ul>"
        "Requires <b>sibc</b> sideset in workpiece .vol + a separate "
        "pre-meshed coil <b>.vol</b> (source/sink sidesets on caps)."
    ),
    METHOD_PEEC_FEM_KELVIN: (
        "<b>PEEC coil + FEM wp + Kelvin</b> (~4-8 min).<br>"
        "<ul>"
        "<li>Coil as PEEC filament (no coil mesh, STEP input)</li>"
        "<li>Workpiece as volumetric FEM with SIBC Robin BC on surface</li>"
        "<li>Kelvin exterior domain for open boundary</li>"
        "<li>Total-field line-integral RHS (scattered mode retired 2026-04-24)</li>"
        "</ul>"
        "Requires <b>sibc</b> sideset + <b>kelvin</b> material in .vol "
        "+ <b>STEP</b> coil.  NO coil material / source / sink needed."
    ),
    METHOD_FEM_FULL: (
        "<b>Full L + P_wp + P_coil</b> (~1-7 min depending on mesh).<br>"
        "<ul>"
        "<li>A-V compound FES (HCurl A + H1 phi on coil)</li>"
        "<li>Coil volumetric + SIBC Robin BC on wp surface</li>"
        "<li>Kelvin exterior for open boundary</li>"
        "<li>L captures wp back-reaction (Lenz/ferromagnetic push)</li>"
        "</ul>"
        "Requires <b>coil</b> material + <b>source</b>/<b>sink</b>/"
        "<b>sibc</b>/<b>kelvin</b> in .vol."
    ),
    METHOD_THERMAL_3D_STATIC: (
        "<b>Thermal: 3D static</b> -- Phase B, no rotation.<br>"
        "<ul>"
        "<li>Reads q_surf .sol from a previous EM solve "
        "(calc_fem_kelvin or other PEEC+FEM methods)</li>"
        "<li>3D volumetric heat equation on a SEPARATE workpiece "
        "thermal mesh (workpiece as a real solid, not a hole)</li>"
        "<li>q_surf held azimuthally fixed; workpiece stationary "
        "(rotation_rpm = 0)</li>"
        "<li>Use for: static one-shot heat-up, feasibility study, "
        "non-rotating IH (induction welding, brazing)</li>"
        "</ul>"
        "Drives calc_heat.py."
    ),
    METHOD_THERMAL_3D_ROTATING: (
        "<b>Thermal: 3D + rotation</b> -- Phase B, true rotation.<br>"
        "<ul>"
        "<li>3D volumetric heat equation + workpiece body spinning "
        "around +z axis</li>"
        "<li>q_surf re-projected on the body frame at each timestep "
        "(v4.58.0+: ``(x*cosθ - y*sinθ, x*sinθ + y*cosθ, z)``); "
        "mesh / FES / mass / stiffness held fixed, only LinearForm "
        "RHS reassembles per step</li>"
        "<li>Per-step overhead ~10 ms on typical meshes</li>"
        "<li>Use for: non-axisymmetric workpiece OR non-axisymmetric "
        "coil with rotation (general spinning case)</li>"
        "</ul>"
        "Drives calc_heat.py with rotation_rpm > 0."
    ),
    METHOD_THERMAL_AXISYM: (
        "<b>Thermal: 2D axisymmetric</b> -- Phase B, axisym shortcut."
        "<br>"
        "<ul>"
        "<li>Rotationally-symmetric workpiece (cylinder, stepped "
        "shaft, disk) meshed in the (r, z) plane</li>"
        "<li>10-100× faster than equivalent 3D mesh</li>"
        "<li>Cross-mesh q_surf transfer is φ-averaged so a "
        "slightly non-axisymmetric coil (gapped torus) still "
        "produces a physically sensible q</li>"
        "<li>rotation_rpm is recorded as metadata (rotation is "
        "implicit in the axisym assumption)</li>"
        "<li>Use for: continuous rotating IH of cylindrical "
        "workpieces (the common Kubota / Kameari workflow)</li>"
        "</ul>"
        "Drives calc_heat_axisym.py."
    ),
}


# ============================================================
# Material presets (common conductors for IH)
# ============================================================

# name -> {sigma [S/m], mu_r, color, default_for_wp, default_for_coil}
MATERIAL_PRESETS = {
    "Copper":          {"sigma": 5.8e7, "mu_r": 1.0},
    "Aluminum":        {"sigma": 3.5e7, "mu_r": 1.0},
    "Brass":           {"sigma": 1.5e7, "mu_r": 1.0},
    "Steel (mu_r=100)": {"sigma": 5.0e6, "mu_r": 100.0},
    "Steel (mu_r=500)": {"sigma": 5.0e6, "mu_r": 500.0},
    "Stainless 304":   {"sigma": 1.4e6, "mu_r": 1.0},
    "Custom":          {"sigma": None,  "mu_r": None},
}


MU_0 = 4e-7 * math.pi


def skin_depth(freq_hz, sigma, mu_r=1.0):
    omega = 2 * math.pi * freq_hz
    return math.sqrt(2.0 / (omega * mu_r * MU_0 * sigma))


# ============================================================
# .vol label inspection (requires NGSolve lazily)
# ============================================================

def inspect_vol_labels(vol_path):
    """Return (materials, boundaries) sets from the .vol.  Returns
    ``(None, None)`` if the file cannot be loaded."""
    if not vol_path or not os.path.isfile(vol_path):
        return None, None
    try:
        # NGSolve import is heavy, do it lazily.
        from ngsolve import Mesh
        m = Mesh(vol_path)
        return set(m.GetMaterials()), set(m.GetBoundaries())
    except Exception:
        return None, None


def _meaningful_labels(labels, *, max_n=12):
    """Return user-named labels (filtering out Cubit auto-names like
    'Surface_42', 'volume_17', 'boundary' so typos like 'sorce' stand
    out at a glance).  Cap at ``max_n`` items to keep the validation
    line short.
    """
    import re
    if not labels:
        return []
    auto_pat = re.compile(r"^(Surface_|volume_|edge_|curve_)\d+$|^boundary$")
    user = sorted({lbl for lbl in labels if not auto_pat.match(lbl)})
    if len(user) > max_n:
        user = user[:max_n] + [f"... (+{len(labels) - max_n} more)"]
    return user


def _label_hint(present, kind):
    """Format a 'available <kind>: [...]' suffix listing the meaningful
    labels actually in the .vol.  Empty string when there are no
    user-named labels (the .vol is generic; nothing useful to list).
    """
    user = _meaningful_labels(present)
    if not user:
        return ""
    return f" available {kind}: {user}"


def check_method_requirements(method, mats, bnds):
    """Return (ok: bool, errors: list[str], warnings: list[str])."""
    errors = []
    warnings = []
    # Vacuum inductance modes (coil only, peec or bem-a) do NOT use .vol.
    if method in (METHOD_PEEC_IND, METHOD_BEMA_IND):
        return True, errors, warnings

    if mats is None:
        # .vol could not be loaded or not yet selected — skip silent.
        return True, errors, warnings

    mat_hint = _label_hint(mats, "materials")
    bnd_hint = _label_hint(bnds, "boundaries")

    # Far-field truncation: NO panel-side warning.
    #
    # The user's intent (Kelvin vs regular FEM with truncation) is
    # captured at .vol export time by the launcher's "Add Kelvin
    # open boundary (auto)" checkbox.  When unchecked the user has
    # explicitly opted into a regular FEM run -- the panel sees the
    # .vol after the fact and cannot tell intent from accident, so
    # warning here just adds noise to a deliberate choice
    # (2026-04-28 user feedback).
    #
    # The audit trail still exists: calc_fem_kelvin emits a
    # ``FARFIELD:`` log line per run telling exactly which path
    # (Periodic Kelvin / outer Dirichlet / gauge-only) actually
    # fired.  Users who care can read the panel debug log; users
    # who don't are not interrupted.
    #
    # Empirical reminder: HCurl(dirichlet="GND_vertex_tag") is a
    # NO-OP because HCurl DOFs live on edges.  Only FACE Dirichlet
    # ("outer") actually constrains HCurl A.  See calc_fem_kelvin
    # FARFIELD log for which case fires.

    if method == METHOD_FEM_FULL:
        if "coil" not in mats:
            errors.append(
                f"Missing material 'coil' (coil volume).{mat_hint}")
        for b in ("source", "sink", "sibc"):
            if b not in bnds:
                errors.append(
                    f"Missing boundary '{b}'.  FEM A-V needs gap-face "
                    f"ports (source/sink) + wp surface (sibc)."
                    f"{bnd_hint}")
    elif method == METHOD_PEEC_FEM_KELVIN:
        # PEEC coil + FEM wp: no coil material, no source/sink, just the
        # workpiece SIBC surface (Kelvin is intent-driven via the
        # launcher checkbox; the panel does not pre-warn here).
        if "sibc" not in bnds:
            errors.append(
                f"Missing boundary 'sibc' (workpiece surface)."
                f"{bnd_hint}")
    else:  # PEEC+BEM
        if "sibc" not in bnds:
            errors.append(
                f"Missing boundary 'sibc' (workpiece surface)."
                f"{bnd_hint}")

    return (len(errors) == 0), errors, warnings


# ============================================================
# IH Panel
# ============================================================

class IHPanel(ModePanel):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._vol_mats = None   # populated by window on .vol load
        self._vol_bnds = None
        # Status messages surfaced by _update_status when a widget
        # value was silently clamped by a method switch or a
        # cross-session restore.  None = no message.
        self._fes_clamp_msg = None
        self._build_ui()

    # ----------------------- UI construction -----------------------

    # _add_section is inherited from ModePanel (hoisted 2026-04-26).

    def _build_ui(self):
        # Method selector
        self._add_section("Method")
        self._method_combo = self.add_combo(
            "method", "Method:",
            [METHOD_PEEC_IND, METHOD_BEMA_IND,
             METHOD_PEEC_BEM, METHOD_BEMA_BEM,
             METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL,
             METHOD_THERMAL_3D_STATIC,
             METHOD_THERMAL_3D_ROTATING,
             METHOD_THERMAL_AXISYM])
        self._method_combo.currentTextChanged.connect(self._on_method_changed)
        self._method_combo.setToolTip(METHOD_TOOLTIP[METHOD_PEEC_IND])

        # Per-method status line (.vol label check / skin depth hint).
        # Base helper: idiom for de-emphasised status row; updates via
        # self._status_label.setText(...) elsewhere in this class.
        self.add_status_label()

        # ============ Drive (frequency + current) ============
        self._add_section("Drive", key="_sec_drive")
        freq = self.add_line("freq", "Frequency [Hz]:", "7000")
        freq.editingFinished.connect(self._update_status)
        freq.setToolTip("Operating frequency.  Skin depth\n"
                         "  delta = sqrt(2 / (w mu sigma))")

        current = self.add_line("current", "Coil current [A, peak]:", "1.0")
        current.setToolTip(
            "Peak (not RMS) coil port current.  Field quantities are "
            "complex phasors; output P is time-averaged (1/2 Re).")

        # ============ Coil material (INDEPENDENT from WP) ============
        self._add_section("Coil material", key="_sec_coil_mat")
        coil_mat = self.add_combo(
            "coil_material", "Preset:",
            list(MATERIAL_PRESETS.keys()), default=0)  # Copper
        coil_mat.currentTextChanged.connect(self._on_coil_material_changed)
        coil_mat.setToolTip(
            "Coil conductor material. Presets set sigma + mu_r. "
            "'Custom' enables manual entry.")

        self.add_line("coil_sigma", "sigma [S/m]:", "5.8e7")

        # ============ PEEC coil input (CAD STEP) ============
        # (Only for PEEC methods, hidden otherwise by _on_method_changed)
        self._add_section("Coil geometry (PEEC)", key="_sec_peec_step")
        step_w = self.add_browse(
            "peec_step", "STEP:",
            filter_str="STEP (*.step *.stp);;All (*)")
        step_w.textChanged.connect(self._emit_validation)
        # "New..." button: write a CoilBuilder racetrack starter .py
        # next to the chosen .step destination AND run it once to
        # materialise the .step the panel needs.  User edits the .py
        # afterwards + re-runs `python coil.py` to refresh the .step.
        self.add_browse_action(
            "peec_step", "New...", self._on_new_coil_template_for_step,
            fixed_width=60)

        # ============ BEM-A coil input (pre-meshed surface .vol) ============
        # BEM-A reads a Cubit / Netgen .vol whose 'source' / 'sink'
        # boundary labels mark the cap faces.  No on-the-fly OCC re-mesh
        # (consistent with CLAUDE.md "Cubit/NGSolve Complete Separation
        # Policy": .vol is the sole computation interface).
        # (Only for BEM-A methods, hidden otherwise by _on_method_changed)
        self._add_section("Coil geometry (BEM-A)", key="_sec_coil_vol")
        coil_vol_w = self.add_browse(
            "coil_vol", "Coil .vol:",
            filter_str="Coil mesh (*.vol);;All (*)")
        coil_vol_w.textChanged.connect(self._emit_validation)
        # BEM-A coil expects sideset (or block) labels for the source/
        # sink port caps.  Defaults match the lab convention "source"
        # and "sink", but some .jou files use abbreviations like
        # "src"/"snk" or capitalised "Source"/"Sink".  Exposing them as
        # editable fields avoids the v4.37.0 "coil source/sink mesh
        # empty" failure for users on a non-default convention.
        self.add_line("coil_source_name", "Source label:",
                       default="source",
                       placeholder="sideset/block name for the +I cap")
        self.add_line("coil_sink_name", "Sink label:",
                       default="sink",
                       placeholder="sideset/block name for the -I cap")

        # ============ Workpiece geometry (.vol) =====================
        # Replaces the legacy top-level Model row in AnalysisWindow
        # (radia 4.35.0+).  Hidden for vacuum-only modes; shown for
        # weak-coupled / FEM modes that need a workpiece mesh.
        self._add_section("Workpiece geometry", key="_sec_wp_vol")
        wp_vol_w = self.add_browse(
            "wp_vol", "Workpiece .vol:",
            filter_str="Netgen Vol (*.vol);;All (*)")
        wp_vol_w.textChanged.connect(self._on_wp_vol_changed_text)

        # ============ Workpiece material (INDEPENDENT from coil) ============
        self._add_section("Workpiece material", key="_sec_wp_material")
        wp_mat = self.add_combo(
            "wp_material", "Preset:",
            list(MATERIAL_PRESETS.keys()), default=3)  # Steel 100
        wp_mat.currentTextChanged.connect(self._on_wp_material_changed)
        wp_mat.setToolTip(
            "Workpiece material. Presets set sigma + mu_r. "
            "'Custom' enables manual entry.")

        self.add_line("wp_sigma", "sigma [S/m]:", "5.0e6")
        self.add_line("mu_r", "mu_r:", "100")
        ht = self.add_line("half_thickness", "half thickness [m]:",
                            "0.0125")
        ht.setToolTip(
            "Half of the workpiece wall thickness for the Dowell "
            "SIBC formula. For solid bulk wp, use min(R_wp, H_wp/2).")

        # ============ Impedance model (Linear SIBC vs ESIM) ============
        self._add_section("Workpiece impedance model", key="_sec_wp_imp")
        imp = self.add_combo(
            "impedance_model", "Model:",
            ["Linear SIBC",
             "Nonlinear ESIM (Karl iteration)"],
            default=0)
        imp.currentTextChanged.connect(self._on_impedance_changed)
        imp.setToolTip(
            "<b>Linear SIBC</b>: Z_s = (1+j) rho/delta (delta includes mu_r; "
            "standard Leontovich, |Z_s| ~ sqrt(mu_r)). "
            "Ok for Cu/Al, and for steel with a constant mu_r.<br>"
            "<b>ESIM</b>: 1D cell problem solves B-H(H) self-consistently "
            "(Karl iteration). Needed when mu_r varies with H (saturated "
            "steel). Requires a BH-file. Wired for PEEC+BEM, "
            "BEM-A+BEM, PEEC+FEM+Kelvin, and Full FEM (since v4.46).")

        # ESIM-only widgets
        self.add_browse("bh_file", "BH file:", default="",
                         filter_str="BH tables (*.txt *.csv);;All (*)")
        self.add_spin("esim_max_iter", "max iter:", 30, 1, 200)
        self.add_line("esim_tol", "tolerance:", "1e-3")
        # Per-DOF + Anderson: IGTE 2026 paper headline (v4.67+/v4.68+).
        # Defaults match the sweep_v2 production setting that the paper
        # reports (per-panel=True, anderson_m=5, relax=0.5).  Default
        # per-panel is OFF for backward compatibility, but kubota-kun
        # workflow flips it on.
        per_panel_check = self.add_check("esim_per_panel",
                        "Per-DOF ESIM (resolves saturation hot-spots; "
                        "IGTE 2026 paper headline)",
                        default=False)
        per_panel_check.setToolTip(
            "When checked: each surface DOF gets its own Z_s computed "
            "from the locally-extracted |H_t|.  Reproduces the IGTE "
            "2026 paper's Fig. 1 heatmap (-22 to -48% disagreement "
            "vs scalar in the BH-knee regime; P_per/P_scalar ~ 0.62 "
            "at I=100 A / f=50 kHz).  Forces Basis order=1 because "
            "calc_inductance's per-DOF |H_t| extractor "
            "(bem_sibc_solver.extract_H_t_per_dof_grad) currently only "
            "supports P1 BIE basis -- selecting p>=2 with per-panel "
            "raises IndexError in v4.67-v4.72.")
        # When per-panel toggles, re-apply visibility hide for FEM-Kelvin /
        # FEM-Full advanced-knob exclusions AND clamp fes_order to 1.
        per_panel_check.stateChanged.connect(
            lambda _: self._on_impedance_changed(
                self.val("impedance_model")))
        anderson_spin = self.add_spin(
            "esim_anderson_m", "Anderson memory m:", 5, 0, 20)
        anderson_spin.setToolTip(
            "Anderson Type-II acceleration history depth.  Production "
            "value m=5 closes the per-DOF dZ noise floor that plain "
            "damped Picard cannot.  Required (not optional) when "
            "Per-DOF ESIM is on AND the workpiece straddles the BH "
            "knee (typical IH surface hardening, I=100-300 A).  "
            "Not used by PEEC+FEM+Kelvin (calc_fem_kelvin.py has no "
            "--esim-anderson-m flag).")
        relax_line = self.add_line("esim_relax", "Karl relax alpha:", "0.5")
        relax_line.setToolTip(
            "Damped Picard relaxation alpha (0..1).  0.5 is the production "
            "default; lower (0.2-0.3) for deep saturation if Anderson "
            "history is also short.  Not used by PEEC+FEM+Kelvin.")

        # ============ Linear solver (method-dependent) ============
        self._add_section("Linear solver", key="_sec_solver")
        solver = self.add_combo("solver", "Solver:", ["pardiso"])
        solver.setToolTip(
            "<b>Inductance / weak-coupled modes</b> "
            "(PEEC-IND, BEM-A-IND, PEEC+BEM, BEM-A+BEM):<br>"
            "&nbsp;&nbsp;Dense LU — workpiece BEM-SIBC assembled densely "
            "(suitable when wp surface DOFs &lt; ~5k)<br>"
            "&nbsp;&nbsp;HACApK — ACA-compressed H-matrix + GMRES "
            "(for wp surface DOFs &gt; ~5k; O(N log N) memory)<br>"
            "<i>Note:</i> For PEEC coil, the perimeter-filament bundle "
            "is ALWAYS dense LU regardless of this selector.  The size "
            "knob gates the WORKPIECE BIE backend (and the BEM-A coil "
            "saddle backend when the coil solver is BEM-A).<br>"
            "<br>"
            "<b>FEM A-V modes</b> (PEEC+FEM+Kelvin, Full FEM):<br>"
            "&nbsp;&nbsp;pardiso — sparse direct (default, fast, memory-heavy)<br>"
            "&nbsp;&nbsp;AMS — Compact AMS+COCR for HCurl p=1 "
            "(low memory; shifted preconditioner internally)<br>"
            "&nbsp;&nbsp;BDDC — preconditioned CG, recommended for p&gt;=2<br>"
            "&nbsp;&nbsp;iccg — generic fallback (Incomplete Cholesky + CG)")

        # ============ Advanced (collapsed by default) ============
        self._add_section("Advanced", key="_sec_advanced")
        n_peri_w = self.add_spin("peec_n_peri",
                                  "PEEC n_peri (perimeter):", 16, 4, 128)
        n_peri_w.setToolTip(
            "Number of filaments placed on the cross-section perimeter "
            "(thin-skin regime).  Typical: 12-24 for circular, 16-32 "
            "for rectangular.  Requires d / skin depth >= 3.")
        self.add_spin("peec_nwinc", "PEEC nwinc (volume grid):", 3, 1, 10)
        self.add_spin("peec_nhinc", "PEEC nhinc (volume grid):", 3, 1, 10)
        # Basis polynomial order.  Meaning depends on method:
        #   PEEC+BEM / BEM-A+BEM -> --h1-order (BEM Lagrange basis on
        #     the workpiece surface; 1 = P1 hat, 2 = Lagrange P2).
        #   PEEC+FEM+Kelvin / Full FEM -> --fes-order (HCurl volume
        #     basis order; 1-3 supported).
        # Same widget for both because the UX ("how smooth should the
        # basis be?") is identical.  P3 is NOT supported on the BEM
        # side (no Lagrange-P3 in-tree assembler) but IS supported on
        # the FEM side -- the calc_inductance.py CLI hard-rejects
        # h1_order=3 with a clear error.
        # NOTE: the GEOMETRY curve order is NOT a panel knob -- it is
        # fixed by the .vol's baked curvedelements (set at Cubit-export
        # time via ``export netgen "f.vol" order N``) and
        # auto-detected by calc_inductance.py from the companion
        # ``.vol.json``.  A post-load ``mesh.Curve(p)`` silently falls
        # back to flat without a CAD callback, so we never expose it.
        self.add_spin("fes_order", "Basis order:", 1, 1, 3)

        # ============ Thermal sub-panel (method=Thermal only) ============
        # HeatPanel from _heat_panel.py embedded as a single sub-widget;
        # all heat-side fields (qsurf .sol, em_vol, wp_vol thermal mesh,
        # material, h_conv / t_ext, time scheme, dt / t_end, probe,
        # rotation_rpm, mesh type 3D vs axisym) live INSIDE this widget
        # and become visible only when method=Thermal.  is_runnable /
        # build_command / wp_vol_path delegate to this sub-panel when
        # method=Thermal (see below).  The standalone radia_heat.py
        # module was removed in radia 4.62.0; this embedded HeatPanel
        # is the sole home for heat analysis.
        self._add_section("Thermal analysis", key="_sec_thermal")
        from radia._heat_panel import HeatPanel
        from PySide6.QtWidgets import QScrollArea, QFrame
        self._heat_panel = HeatPanel(parent=self)
        # Wrap HeatPanel in a QScrollArea with a bounded max height so
        # the parent window stays inside the 2K-monitor budget
        # (panel_qa MAX_HEIGHT_RED = 1700 px).  The heat panel's
        # natural sizeHint is ~1100-1200 px which would push the
        # window past 1800 px; the scroll viewport caps the row at
        # 700 px and adds a vertical scrollbar when needed.  Width
        # is unconstrained since labels need to fit.
        scroll = QScrollArea()
        scroll.setWidget(self._heat_panel)
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(700)
        scroll.setFrameShape(QFrame.NoFrame)
        self._heat_panel_scroll = scroll
        # Embed the scroll-wrapped heat panel as a single full-width
        # row in IHPanel's form layout.
        self._form.addRow(scroll)
        # Track the heat-panel row so _set_row_visible() can collapse
        # it alongside the section header when method != Thermal.
        self._row_indices["_heat_panel_row"] = self._form.rowCount() - 1
        # Forward the sub-panel's validation signal so the parent Run
        # button enables/disables in sync with the sub-panel's
        # is_runnable() state.
        self._heat_panel.validationChanged = self._emit_validation

        # ============ Initial state ============
        # Apply material presets to refresh coil_sigma / wp_sigma / mu_r
        self._on_coil_material_changed(coil_mat.currentText())
        self._on_wp_material_changed(wp_mat.currentText())

        self._method_combo.setCurrentText(METHOD_PEEC_IND)
        self._on_method_changed(METHOD_PEEC_IND)
        self._update_status()

    def _emit_validation(self, *_):
        """Forward edits to AnalysisWindow so Run enables/disables."""
        cb = getattr(self, "validationChanged", None)
        if callable(cb):
            cb()

    # --------------- CoilBuilder template wizard (STEP) ---------------

    def _on_new_coil_template_for_step(self, line_edit):
        """IH-side companion to the EM "New..." wizard.

        Writes a CoilBuilder starter .py at <basename>.py and immediately
        runs it once to materialise <basename>.step (which the IH PEEC
        modes consume).  Sets the peec_step field to the .step path.

        For later customisation: the user edits the sibling .py and
        re-runs `python <basename>.py` to refresh the .step (the
        template's __main__ writes the .step automatically).
        """
        from PySide6.QtWidgets import QFileDialog, QMessageBox
        from radia_gui_base import COIL_TEMPLATE
        existing = line_edit.text().strip()
        if existing and os.path.isdir(os.path.dirname(existing)):
            suggested = os.path.join(
                os.path.dirname(existing), "coil_new.step")
        else:
            suggested = os.path.abspath("coil.step")
        step_path, _ = QFileDialog.getSaveFileName(
            self, "Save new coil (.py + .step)", suggested,
            "STEP (*.step *.stp);;All (*)")
        if not step_path:
            return
        # Sibling .py path for editing.  COIL_TEMPLATE's __main__
        # writes the .step using `splitext(__file__)[0] + '.step'`,
        # so the .py basename must match the .step basename.
        py_path = os.path.splitext(step_path)[0] + ".py"
        try:
            with open(py_path, "w", encoding="utf-8") as f:
                f.write(COIL_TEMPLATE)
        except OSError as exc:
            QMessageBox.warning(
                self, "Could not write coil .py",
                f"Failed to write {py_path}:\n{exc}")
            return
        # Run the .py in-process so the .step is materialised
        # immediately.  __main__ guards the write_step call so any
        # netgen.occ failure surfaces here without leaving a stale
        # peec_step pointing at a missing file.
        try:
            import runpy
            runpy.run_path(py_path, run_name="__main__")
        except Exception as exc:
            QMessageBox.warning(
                self, "Could not run coil .py",
                f"Wrote {py_path} but failed to materialise .step:\n"
                f"{exc}\n\nEdit the .py and run "
                f"`python {os.path.basename(py_path)}` manually.")
            return
        if not os.path.isfile(step_path):
            QMessageBox.warning(
                self, "STEP not produced",
                f"{py_path} ran but {step_path} was not written.\n"
                "STEP export needs netgen.occ; check the .py output.")
            return
        line_edit.setText(step_path)

    # ----------------------- Material presets -----------------------

    def _apply_material(self, preset_name, sigma_key, mu_r_key):
        preset = MATERIAL_PRESETS.get(preset_name, {})
        sigma = preset.get("sigma")
        mu_r = preset.get("mu_r")
        sigma_w = self._widgets.get(sigma_key)
        mu_w = self._widgets.get(mu_r_key)
        is_custom = (preset_name == "Custom")
        if sigma_w is not None:
            if sigma is not None:
                sigma_w.setText(f"{sigma:.3g}")
            sigma_w.setEnabled(is_custom)
        if mu_w is not None:
            if mu_r is not None:
                mu_w.setText(f"{mu_r:.3g}")
            mu_w.setEnabled(is_custom)

    def _on_coil_material_changed(self, name):
        self._apply_material(name, "coil_sigma", None)
        self._update_status()

    def _on_wp_material_changed(self, name):
        self._apply_material(name, "wp_sigma", "mu_r")
        # Hint: Steel materials often need ESIM for accurate saturation
        # prediction.  Don't auto-switch — just nudge via status line.
        self._update_status()

    def _on_impedance_changed(self, name):
        is_esim = name.startswith("Nonlinear ESIM")
        for key in ("bh_file", "esim_max_iter", "esim_tol",
                    "esim_per_panel", "esim_anderson_m", "esim_relax"):
            self._set_row_visible(key, is_esim)
        # Re-apply the PEEC_FEM_KELVIN esim_tol hide (it's the one
        # ESIM mode that has no --esim-tol equivalent on the calc
        # side, so we must keep that row invisible even when ESIM is
        # selected).  Same applies to anderson-m and relax which
        # calc_fem_kelvin.py also does not accept (verified
        # 2026-05-24 by panel-cli-diff against calc_fem_kelvin.py).
        if is_esim and self.val("method") == METHOD_PEEC_FEM_KELVIN:
            self._set_row_visible("esim_tol", False)
            self._set_row_visible("esim_anderson_m", False)
            self._set_row_visible("esim_relax", False)
        # Full-FEM (calc_fem_coilmesh.py) accepts per-panel + relax
        # but NOT anderson-m (verified same audit).
        if is_esim and self.val("method") == METHOD_FEM_FULL:
            self._set_row_visible("esim_anderson_m", False)
        # Per-DOF ESIM only works at h1-order=1 in v4.67-v4.72:
        # extract_H_t_per_dof_grad in bem_sibc_solver.py assumes phi_vec
        # is the P1 vertex-DOF vector and IndexErrors at p>=2.  Clamp
        # fes_order max to 1 when per-panel is checked.  When unchecked,
        # restore the method-dependent max directly (mirrors the logic
        # in _on_method_changed lines 875-878 -- do NOT call
        # _on_method_changed recursively because it would re-fire
        # _on_impedance_changed and loop).
        fes_spin = self._widgets.get("fes_order")
        per_on = (is_esim and self._widgets["esim_per_panel"].isChecked())
        if fes_spin is not None:
            if per_on and fes_spin.maximum() != 1:
                fes_spin.setMaximum(1)
                self._fes_clamp_msg = (
                    "fes_order clamped to 1 because Per-DOF ESIM "
                    "requires P1 BIE basis (extract_H_t_per_dof_grad "
                    "P2+ support is open).")
            elif not per_on:
                method = self._method_combo.currentText()
                is_weak_now = method in (METHOD_PEEC_BEM, METHOD_BEMA_BEM)
                is_fem_now = method in (METHOD_PEEC_FEM_KELVIN,
                                          METHOD_FEM_FULL)
                if is_weak_now or is_fem_now:
                    fes_spin.setMaximum(2 if is_weak_now else 3)
        self._update_status()

    def restore_state(self, state):
        """Override to detect cross-session fes_order silent clamp.

        ModePanel.restore_state walks widgets in dict-insertion order;
        the method combo is restored before fes_order and its
        currentTextChanged signal fires _on_method_changed which sets
        the fes_order spin's max to 2 for weak-coupled methods.  When
        fes_order=3 was saved against a weak method, Qt then silently
        clamps the value to 2.  We compare the saved value to the
        actual restored value and surface a status warning when they
        differ.
        """
        super().restore_state(state)
        if not state:
            return
        saved_fes = state.get("fes_order")
        if saved_fes is None or "fes_order" not in self._widgets:
            return
        try:
            target = int(saved_fes)
        except (ValueError, TypeError):
            return
        actual = self._widgets["fes_order"].value()
        if target > actual:
            # Append to whatever _on_method_changed already set.
            note = (f"fes_order saved={target} -> restored={actual} "
                    f"(current method max).  Switch to a non-weak "
                    f"method to use p={target} again.")
            self._fes_clamp_msg = (
                note if not self._fes_clamp_msg
                else f"{self._fes_clamp_msg}; {note}")
            self._update_status()

    # ----------------------- Method + visibility -----------------------

    def _on_method_changed(self, method):
        is_peec_ind = (method == METHOD_PEEC_IND)
        is_bema_ind = (method == METHOD_BEMA_IND)
        is_peec_bem = (method == METHOD_PEEC_BEM)
        is_bema_bem = (method == METHOD_BEMA_BEM)
        is_peec_fem_k = (method == METHOD_PEEC_FEM_KELVIN)
        is_fem = (method == METHOD_FEM_FULL)
        is_thermal_3d_static   = (method == METHOD_THERMAL_3D_STATIC)
        is_thermal_3d_rotating = (method == METHOD_THERMAL_3D_ROTATING)
        is_thermal_axisym      = (method == METHOD_THERMAL_AXISYM)
        is_thermal = (is_thermal_3d_static
                       or is_thermal_3d_rotating
                       or is_thermal_axisym)
        # Inductance-style (vacuum) vs weak-coupled (workpiece) groupings.
        is_vacuum = is_peec_ind or is_bema_ind
        is_weak = is_peec_bem or is_bema_bem
        is_inductance_path = is_vacuum or is_weak    # any calc_inductance.py mode
        is_bem_a = is_bema_ind or is_bema_bem
        is_peec = is_peec_ind or is_peec_bem
        # Thermal mode delegates ALL field visibility to the embedded
        # HeatPanel sub-widget -- the EM-side sections (drive, coil
        # material, coil geometry, workpiece material, ...) are not
        # meaningful here and stay hidden.  When method != Thermal the
        # heat panel section hides.
        is_em = not is_thermal
        # PEEC + PEEC+FEM+Kelvin take CAD STEP (filament coil); BEM-A
        # variants take a pre-meshed .vol coil.  FEM-full uses a
        # volumetric coil baked into its workpiece .vol.
        needs_step = is_em and (is_peec or is_peec_fem_k)
        needs_coil_vol = is_em and is_bem_a
        needs_wp = is_em and (is_weak or is_peec_fem_k or is_fem)

        self._method_combo.setToolTip(METHOD_TOOLTIP.get(method, ""))

        # Solver selection per method.  Always visible — users want to
        # know what is solving their problem.
        solver = self._widgets["solver"]
        prev = solver.currentText()
        solver.clear()
        if is_inductance_path:
            # PEEC bundle solver (peec) or BEM-A coil-saddle solver (bem-a).
            # Both expose the same Dense/HACApK choice for the workpiece BIE.
            solver.addItems(["Dense LU (small)",
                              "HACApK (large)"])
        else:  # FEM side: PEEC+FEM+Kelvin and FEM-full share solver choices
            solver.addItems(["pardiso (direct)",
                              "AMS (iterative, p=1)",
                              "BDDC (iterative, p>=2)",
                              "iccg (fallback)"])
        idx = solver.findText(prev)
        if idx >= 0:
            solver.setCurrentIndex(idx)

        # Visibility per method (rows + their section headers):
        # PEEC modes -> CAD STEP row; BEM-A modes -> pre-meshed .vol row.
        self._set_row_visible("_sec_peec_step", needs_step)
        self._set_row_visible("peec_step", needs_step)
        self._set_row_visible("_sec_coil_vol", needs_coil_vol)
        self._set_row_visible("coil_vol", needs_coil_vol)
        # Source/sink labels only meaningful for BEM-A coil .vol.
        self._set_row_visible("coil_source_name", needs_coil_vol)
        self._set_row_visible("coil_sink_name", needs_coil_vol)
        # Perimeter placement (n_peri) only for PEEC variants.  BEM-A
        # variants no longer have a runtime mesh-maxh knob (the .vol
        # was meshed externally by Cubit / Netgen).
        self._set_row_visible("peec_n_peri", is_peec)
        # Volume filament grid (nwinc/nhinc) only PEEC+FEM+Kelvin now.
        self._set_row_visible("peec_nwinc", is_peec_fem_k)
        self._set_row_visible("peec_nhinc", is_peec_fem_k)
        # fes_order spin is shown for all weak-coupled / FEM modes:
        #   * PEEC+BEM / BEM-A+BEM    -> --h1-order (Lagrange basis order)
        #   * PEEC+FEM+Kelvin / Full FEM -> --fes-order (HCurl volume order)
        # Hidden for vacuum-only modes (no workpiece, no FEM volume).
        # Note: for BEM modes, this is the BASIS order only.  The
        # geometry CURVE order is fixed at Cubit-export time
        # (.vol.json) and auto-detected by calc_inductance.py -- a
        # post-load mesh.Curve(p) silently falls back to flat, so we
        # don't expose it as a knob.
        self._set_row_visible("fes_order",
                              is_weak or is_peec_fem_k or is_fem)
        # Weak-coupled BEM modes flow through calc_inductance.py
        # --h1-order which has choices=[1, 2] -- selecting 3 would
        # fire an argparse error.  FEM-side modes (PEEC+FEM+Kelvin,
        # FEM-full) accept --fes-order up to 3.
        fes_spin = self._widgets["fes_order"]
        new_max = 2 if is_weak else 3
        old_val = fes_spin.value()
        fes_spin.setMaximum(new_max)
        # Qt silently clamps spin.value() when setMaximum drops it
        # below the current value.  Surface the demotion in the
        # status line so the user sees it (otherwise switching from
        # FEM-full p=3 to PEEC_BEM would silently demote to p=2).
        if old_val > new_max:
            self._fes_clamp_msg = (
                f"fes_order {old_val} -> {new_max} for {method} "
                f"(weak-coupled BEM is choices=[1,2])")
        else:
            self._fes_clamp_msg = None
        # Workpiece geometry (.vol) is shown for all modes that need a
        # workpiece mesh; the legacy top-level Model row is gone.
        self._set_row_visible("_sec_wp_vol", needs_wp)
        self._set_row_visible("wp_vol", needs_wp)
        self._set_row_visible("_sec_wp_material", needs_wp)
        self._set_row_visible("_sec_wp_imp", needs_wp)

        # Workpiece widgets are meaningless for vacuum-only modes.
        for key in ("wp_material", "wp_sigma", "mu_r", "half_thickness",
                    "impedance_model", "bh_file",
                    "esim_max_iter", "esim_tol",
                    "esim_per_panel", "esim_anderson_m", "esim_relax"):
            self._set_row_visible(key, needs_wp)

        # ESIM sub-widgets re-evaluated by _on_impedance_changed when
        # needs_wp is True.
        if needs_wp:
            self._on_impedance_changed(self.val("impedance_model"))

        # esim_tol is forwarded via --esim-tol by calc_inductance.py
        # (PEEC_BEM / BEMA_BEM) and calc_fem_coilmesh.py (FEM_FULL),
        # but calc_fem_kelvin.py (PEEC_FEM_KELVIN) only accepts
        # --max-iter -- it has no --esim-tol equivalent.  Hiding the
        # row in PEEC_FEM_KELVIN avoids the user setting a tolerance
        # that the solver never sees.
        if is_peec_fem_k:
            self._set_row_visible("esim_tol", False)

        # Thermal mode: hide EVERY EM-side section + row and show the
        # embedded HeatPanel.  When switching back to an EM method
        # those sections re-enable themselves through the per-method
        # branches above (needs_step / needs_wp / ...).
        for em_sec in ("_sec_drive", "_sec_coil_mat", "_sec_solver",
                        "_sec_advanced"):
            self._set_row_visible(em_sec, not is_thermal)
        for em_row in ("freq", "current", "coil_material", "coil_sigma",
                        "solver"):
            self._set_row_visible(em_row, not is_thermal)
        self._set_row_visible("_sec_thermal", is_thermal)
        self._set_row_visible("_heat_panel_row", is_thermal)
        # The HeatPanel's own widgets sit inside a sub-form; toggling
        # the parent row hides the whole sub-panel.  Direct setVisible
        # is the safety belt for Qt versions where setRowVisible only
        # collapses the QFormLayout row, not the embedded widget.
        # The scroll wrapper also needs toggling so its viewport
        # doesn't contribute to sizeHint when thermal mode is off.
        if hasattr(self, "_heat_panel_scroll"):
            self._heat_panel_scroll.setVisible(is_thermal)
        if hasattr(self, "_heat_panel"):
            self._heat_panel.setVisible(is_thermal)
            # Auto-set the embedded HeatPanel's mesh_type +
            # rotation_rpm visibility based on which Thermal method
            # is active.  The mesh_type combo is HIDDEN (the parent
            # Method dropdown encodes that choice); rotation_rpm is
            # also hidden + reset to 0 for the 3D-static method.
            from radia._heat_panel import (
                MESH_TYPE_3D, MESH_TYPE_AXISYM,
            )
            mesh_w = self._heat_panel._widgets.get("mesh_type")
            rpm_w = self._heat_panel._widgets.get("rotation_rpm")
            n_phi_w = self._heat_panel._widgets.get("n_phi_samples")
            if is_thermal_axisym and mesh_w is not None:
                mesh_w.setCurrentText(MESH_TYPE_AXISYM)
            elif (is_thermal_3d_static or is_thermal_3d_rotating) \
                    and mesh_w is not None:
                mesh_w.setCurrentText(MESH_TYPE_3D)
            # Hide the mesh_type row -- method dropdown owns that choice.
            # Also hide the "Mesh type" section header that sits above it,
            # otherwise it becomes an orphan (no content under it).
            self._heat_panel._set_row_visible("_sec_mesh_type", False)
            self._heat_panel._set_row_visible("mesh_type", False)
            # rotation_rpm: hidden + zeroed for static, visible for
            # rotating + axisym.  Cache the user-entered value before
            # zeroing so a 3D-rotating(1200) -> 3D-static -> 3D-rotating
            # round-trip restores 1200 instead of leaving the user to
            # retype.  Bug 2026-05-25: the previous code unconditionally
            # wrote "0" into rpm_w on every static-method visit and lost
            # whatever the user had typed before.
            #
            # rotation_rpm is a QDoubleSpinBox since v4.78.0 (was
            # QLineEdit pre-v4.78.0); use .value() / .setValue() on it.
            from PySide6.QtWidgets import QDoubleSpinBox, QLineEdit
            def _rpm_get():
                if isinstance(rpm_w, QDoubleSpinBox):
                    return float(rpm_w.value())
                if isinstance(rpm_w, QLineEdit):
                    txt = rpm_w.text().strip()
                    return float(txt) if txt else 0.0
                return 0.0
            def _rpm_set(v):
                if isinstance(rpm_w, QDoubleSpinBox):
                    rpm_w.setValue(float(v))
                elif isinstance(rpm_w, QLineEdit):
                    rpm_w.setText(str(v))
            if rpm_w is not None:
                if is_thermal_3d_static:
                    current = _rpm_get()
                    if current > 0.0:
                        # User had a non-zero rpm before switching to
                        # static -- remember it so we can restore on
                        # return to a rotating-capable mode.
                        self._heat_panel.setProperty(
                            "_last_rotation_rpm", current)
                    _rpm_set(0.0)
                    self._heat_panel._set_row_visible(
                        "rotation_rpm", False)
                else:
                    self._heat_panel._set_row_visible(
                        "rotation_rpm", is_thermal_3d_rotating
                                         or is_thermal_axisym)
                    # Restore cached value when returning from static
                    # (only if the user hasn't already typed something).
                    if _rpm_get() == 0.0:
                        cached = self._heat_panel.property(
                            "_last_rotation_rpm")
                        if cached is not None:
                            try:
                                _rpm_set(float(cached))
                            except (TypeError, ValueError):
                                pass
            # axisym: rotation_rpm is metadata-only -- the answer does
            # NOT depend on rpm in the axisym solver (rotation is
            # implicit in the axisymmetric assumption).  Grey out the
            # widget so the user sees that typing a value won't change
            # the result.  Also grey rotation_axis (rotation is around
            # the implicit symmetry axis, which the .vol's r/z plane
            # encodes -- no axis choice to make).
            axis_w = self._heat_panel._widgets.get("rotation_axis")
            if is_thermal_axisym:
                if rpm_w is not None:
                    rpm_w.setEnabled(False)
                    rpm_w.setToolTip(
                        "Axisym: rotation is implicit (workpiece is "
                        "rotation-symmetric by construction); the "
                        "value is recorded as metadata only and does "
                        "NOT affect the solve.  Pick a 3D Thermal "
                        "method if you want rpm to drive the physics.")
                if axis_w is not None:
                    axis_w.setEnabled(False)
                    self._heat_panel._set_row_visible(
                        "rotation_axis", False)
            else:
                # Non-axisym thermal methods restore the rotation_axis
                # row + re-enable both widgets (unless Source=Uniform
                # has already disabled them via _on_heat_source_changed).
                if axis_w is not None:
                    self._heat_panel._set_row_visible(
                        "rotation_axis",
                        is_thermal_3d_rotating)
                # Re-trigger source change to refresh Uniform/Spatial
                # enable state (Uniform also greys these widgets).
                src_w = self._heat_panel._widgets.get("heat_source")
                if src_w is not None:
                    self._heat_panel._on_heat_source_changed(
                        src_w.currentText())
            # n_phi_samples (axisym only) -- HeatPanel's own
            # _on_mesh_type_changed handler already toggles this
            # based on mesh_type; nothing extra needed here.

        self._update_status()
        self._emit_validation()

    # ----------------------- Status line -----------------------

    def _update_status(self):
        """Compose a short status string: .vol label check + skin depth."""
        lines = []
        method = self._method_combo.currentText()
        ok, errors, warnings = check_method_requirements(
            method, self._vol_mats, self._vol_bnds)

        if self._vol_mats is not None:
            if ok and not warnings:
                lines.append("<span style='color:#080;'>.vol OK</span>")
            if warnings:
                for w in warnings:
                    lines.append(f"<span style='color:#A80;'>warn: {w}</span>")
            if errors:
                for e in errors:
                    lines.append(f"<span style='color:#C00;'>ERROR: {e}</span>")

        # Physics sanity: delta_wp vs half_thickness (workpiece only).
        if method != METHOD_PEEC_IND:
            try:
                freq = float(self.val("freq"))
                wp_sigma = float(self.val("wp_sigma"))
                mu_r = float(self.val("mu_r"))
                half_thickness = float(self.val("half_thickness"))
                delta = skin_depth(freq, wp_sigma, mu_r)
                ratio = half_thickness / delta
                lines.append(
                    f"WP delta = {delta*1e3:.3f} mm, R/delta = {ratio:.1f}")
                if ratio < 3:
                    lines.append(
                        "<span style='color:#A80;'>warn: R/delta &lt; 3; "
                        "SIBC may under-estimate (volumetric FEM preferred)."
                        "</span>")
            except (ValueError, ZeroDivisionError):
                pass

        # Surface any spin-clamp warning (cross-session restore or
        # in-session method switch demoted fes_order silently).
        if self._fes_clamp_msg:
            lines.append(
                f"<span style='color:#A80;'>warn: "
                f"{self._fes_clamp_msg}</span>")

        self._status_label.setText("<br>".join(lines))

    # ----------------------- External hooks -----------------------

    def is_runnable(self):
        method = self._method_combo.currentText()
        # Thermal mode delegates to the embedded HeatPanel which has
        # its own is_runnable() (wp_vol + qsurf .sol + em_vol checks).
        if method in THERMAL_METHODS:
            return self._heat_panel.is_runnable() \
                if hasattr(self, "_heat_panel") else False
        ok, _errors, _warnings = check_method_requirements(
            method, self._vol_mats, self._vol_bnds)
        # Vacuum inductance modes need a coil input only; the workpiece
        # .vol row is hidden.  PEEC takes STEP; BEM-A takes pre-meshed .vol.
        if method == METHOD_PEEC_IND:
            step = self.val("peec_step") if "peec_step" in self._widgets else ""
            return ok and bool(step) and os.path.isfile(step)
        if method == METHOD_BEMA_IND:
            cv = self.val("coil_vol") if "coil_vol" in self._widgets else ""
            return ok and bool(cv) and os.path.isfile(cv)
        return ok and (self._vol_mats is not None)

    # AnalysisWindow._on_run() asks the panel for its workpiece .vol so it
    # can pass an absolute path into build_command and set the subprocess
    # work_dir.  Read from the per-panel wp_vol widget (4.35.0+).  Empty
    # string for vacuum modes that don't have a workpiece.
    def wp_vol_path(self):
        # Thermal mode reads the wp .vol from the embedded HeatPanel.
        if self._method_combo.currentText() in THERMAL_METHODS \
                and hasattr(self, "_heat_panel"):
            return self._heat_panel.wp_vol_path()
        if "wp_vol" not in self._widgets:
            return ""
        return self.val("wp_vol")

    def _on_wp_vol_changed_text(self, _text):
        """Re-inspect the workpiece .vol labels after the user edits the
        wp_vol field (Browse... or manual edit).  Mirrors the legacy
        window-level .vol hook against the panel-owned widget.
        """
        path = self.val("wp_vol")
        # Reset label cache on empty path (e.g. user cleared the field).
        if not path or not os.path.isfile(path):
            self._vol_mats = None
            self._vol_bnds = None
        else:
            mats, bnds = inspect_vol_labels(path)
            self._vol_mats = mats
            self._vol_bnds = bnds
        self._update_status()
        self._emit_validation()

    # ----------------------- Command building -----------------------

    def build_command(self, vol_path):
        method = self.val("method")
        # Thermal mode: delegate to the embedded HeatPanel's
        # build_command (calc_heat.py / calc_heat_axisym.py based on
        # mesh_type which the parent _on_method_changed has already
        # set from the Method dropdown choice).  vol_path here is
        # the HeatPanel's own --wp-vol; build_command on HeatPanel
        # ignores it and reads its own ``wp_vol`` widget.
        if method in THERMAL_METHODS:
            return self._heat_panel.build_command(vol_path)
        # Vacuum modes (no .vol).  Same builder for PEEC and BEM-A;
        # _coil_solver_cli() picks the --coil-solver flag from the
        # method name.
        if method in (METHOD_PEEC_IND, METHOD_BEMA_IND):
            return self._build_peec_inductance_command()
        if not vol_path:
            raise ValueError("No .vol file specified.")
        # Weak-coupled modes (workpiece via scalar BEM-SIBC).
        if method in (METHOD_PEEC_BEM, METHOD_BEMA_BEM):
            return self._build_peec_bem_command(vol_path)
        if method == METHOD_PEEC_FEM_KELVIN:
            return self._build_fem_kelvin_command(vol_path)
        return self._build_fem_coilmesh_command(vol_path)

    # UI solver text -> CLI arg mapping
    # Inductance / weak-coupled modes: one "size" knob maps to both the
    # coil-side BEM solver (--coil-bem-solver) and the workpiece-side
    # BEM backend (--wp-bem-backend).  Dense LU is the small/exact path;
    # HACApK is the large/ACA-compressed path.
    _PEEC_SOLVER_MAP = {
        "Dense LU (small)": {"coil_bem_solver": "dense-lu",
                              "wp_bem_backend":  "intree-dense"},
        "HACApK (large)":   {"coil_bem_solver": "hacapk-gmres",
                              "wp_bem_backend":  "hacapk"},
    }
    _FEM_SOLVER_MAP = {
        "pardiso (direct)":              "pardiso",
        "AMS (iterative, p=1)":          "ams",
        "BDDC (iterative, p>=2)":        "bddc",
        "iccg (fallback)":               "iccg",
    }

    def _impedance_model_cli(self):
        imp = self.val("impedance_model")
        return "esim" if imp.startswith("Nonlinear ESIM") else "sibc"

    def _coil_solver_cli(self):
        """Derive the ``--coil-solver`` flag from the active method.

        BEM-A variants of inductance / weak-coupled modes use the
        BEM-A coil saddle solver; PEEC variants use the PEEC perimeter
        bundle solver.  Other methods (PEEC+FEM+Kelvin, FEM-full) do
        not flow through ``calc_inductance.py``.
        """
        method = self._method_combo.currentText() \
            if hasattr(self, "_method_combo") else METHOD_PEEC_IND
        if method in (METHOD_BEMA_IND, METHOD_BEMA_BEM):
            return "bem-a"
        return "peec"

    # CLI-DIFF: ignore --output -- auto-injected by calc_main wrapper
    # (calc_common.py:1173-1177).  Static scanners flag this as REJECT
    # for every builder, but the flag IS accepted at runtime.
    # CLI-DIFF: ignore --bh-file --esim-anderson-m --esim-max-iter --esim-per-panel --esim-relax --esim-tol --h1-order --half-thickness --impedance-model --mu-r --sigma --vol --wp-label --wp-bem-backend -- coil-only mode shares calc_inductance.py but intentionally omits workpiece and ESIM flags.
    # CLI-DIFF: ignore --coil-aca-eps --coil-gmres-tol --coil-maxh --coil-msh-output --coil-only --coil-rwg-quad-degree --coil-rwg-singular-nq --coil-saddle-solver --coupling-mode --n-threads --peec-proximity --peec-proximity-max-iter --peec-proximity-relax --peec-proximity-tol --telegen-form --wp-aca-eps --wp-gmres-tol --write-summary -- expert/diagnostic calc_inductance.py knobs are deliberately CLI-only in the production IH panel.

    def _build_peec_inductance_command(self):
        """Coil-only inductance (vacuum, no workpiece).

        Replaces the old ``calc_peec_inductance.py`` invocation by
        calling the unified ``calc_inductance.py`` with ``--coil-solver``
        chosen by the panel.  Method retains the old name for backward
        compatibility with the IH panel registry / hooks; the underlying
        CLI is the new unified one.

        PEEC -> --coil-step (CAD); BEM-A -> --coil-vol (pre-meshed).
        """
        coil_solver = self._coil_solver_cli()
        if coil_solver == "peec":
            coil_in = self.val("peec_step")
            if not coil_in:
                raise ValueError("Coil STEP file is required for PEEC mode.")
            if not os.path.isfile(coil_in):
                raise ValueError(f"STEP file not found: {coil_in}")
            coil_arg = ["--coil-step", coil_in]
        else:  # bem-a
            coil_in = self.val("coil_vol") if "coil_vol" in self._widgets else ""
            if not coil_in:
                raise ValueError("Coil .vol file is required for BEM-A mode.")
            if not os.path.isfile(coil_in):
                raise ValueError(f"Coil .vol not found: {coil_in}")
            coil_arg = ["--coil-vol", coil_in,
                        "--coil-source-name", self.val("coil_source_name"),
                        "--coil-sink-name",   self.val("coil_sink_name")]
        bem_size = self._PEEC_SOLVER_MAP.get(
            self.val("solver"),
            {"coil_bem_solver": "auto", "wp_bem_backend": "hacapk"})
        cmd = [_PYTHON, calc_script("calc_inductance.py"),
               *coil_arg,
               "--coil-solver", coil_solver,
               "--coil-bem-solver", bem_size["coil_bem_solver"],
               "--frequency", self.val("freq"),
               "--current", self.val("current"),
               "--coil-sigma", self.val("coil_sigma"),
               "--msh-output", msh_output(coil_in, "_peec_ind"),
               "--output", json_output(coil_in, "_peec_ind")]
        if coil_solver == "peec":
            cmd += ["--peec-n-peri", str(self.val("peec_n_peri"))]
        return cmd

    def _build_peec_bem_command(self, vol_path):
        """Weak-coupled coil + workpiece (Telegen ΔL via scalar BEM-SIBC).

        Unified CLI: ``calc_inductance.py`` with ``--coil-solver`` +
        ``--vol`` for weak coupling.  Replaces ``calc_peec_bem.py``
        (PEEC coil) and ``calc_coil_bem_a_workpiece.py`` (BEM-A coil)
        in a single dispatch.

        PEEC -> --coil-step (CAD); BEM-A -> --coil-vol (pre-meshed).
        """
        coil_solver = self._coil_solver_cli()
        if coil_solver == "peec":
            coil_in = self.val("peec_step")
            if not coil_in:
                raise ValueError("Coil STEP file is required for PEEC weak-coupled mode.")
            if not os.path.isfile(coil_in):
                raise ValueError(f"STEP file not found: {coil_in}")
            coil_arg = ["--coil-step", coil_in]
        else:  # bem-a
            coil_in = self.val("coil_vol") if "coil_vol" in self._widgets else ""
            if not coil_in:
                raise ValueError("Coil .vol file is required for BEM-A weak-coupled mode.")
            if not os.path.isfile(coil_in):
                raise ValueError(f"Coil .vol not found: {coil_in}")
            coil_arg = ["--coil-vol", coil_in,
                        "--coil-source-name", self.val("coil_source_name"),
                        "--coil-sink-name",   self.val("coil_sink_name")]
        bem_size = self._PEEC_SOLVER_MAP.get(
            self.val("solver"),
            {"coil_bem_solver": "auto", "wp_bem_backend": "hacapk"})
        cmd = [_PYTHON, calc_script("calc_inductance.py"),
               *coil_arg,
               "--coil-solver", coil_solver,
               "--coil-bem-solver", bem_size["coil_bem_solver"],
               "--wp-bem-backend",  bem_size["wp_bem_backend"],
               "--frequency", self.val("freq"),
               "--current", self.val("current"),
               "--coil-sigma", self.val("coil_sigma"),
               "--vol", vol_path,
               "--wp-label", "sibc",
               "--sigma", self.val("wp_sigma"),
               "--half-thickness", self.val("half_thickness"),
               "--mu-r", self.val("mu_r"),
               "--impedance-model", self._impedance_model_cli(),
               "--h1-order", str(self.val("fes_order")),
               "--msh-output", msh_output(vol_path, "_peec_bem"),
               "--output", json_output(vol_path, "_peec_bem")]
        # Note: --h1-order is the BEM Lagrange BASIS order (1 or 2),
        # user-selectable via the "Basis order" spin.  The geometry
        # CURVE order is independent: it is fixed at Cubit-export time
        # (companion .vol.json's "order" field) and auto-detected by
        # calc_inductance.py -- there is no knob for it because a
        # post-load mesh.Curve(p) silently falls back to flat without
        # a CAD callback.
        if coil_solver == "peec":
            cmd += ["--peec-n-peri", str(self.val("peec_n_peri"))]
        if self._impedance_model_cli() == "esim":
            bh = self.val("bh_file")
            if not bh:
                raise ValueError(
                    "Nonlinear ESIM impedance model requires a BH "
                    "file (2-column H[A/m] B[T]).  Fill the 'BH file' "
                    "field, or switch to Linear SIBC.")
            cmd += ["--bh-file", bh,
                    "--esim-max-iter", str(self.val("esim_max_iter")),
                    "--esim-tol", self.val("esim_tol"),
                    "--esim-anderson-m", str(self.val("esim_anderson_m")),
                    "--esim-relax", self.val("esim_relax")]
            if self._widgets["esim_per_panel"].isChecked():
                cmd.append("--esim-per-panel")
        return cmd

    def _build_fem_kelvin_command(self, vol_path):
        """PEEC coil (filament) + FEM workpiece (SIBC Robin) + Kelvin.

        Uses ``calc_fem_kelvin.py --formulation total`` (the only
        surviving formulation after scattered was retired 2026-04-24
        for unfixable P_wp under-prediction).  The coil is driven by
        the PEEC filament line-integral RHS, not by a coil mesh.

        Workpiece material is passed via ``--material custom`` +
        explicit ``--sigma`` / ``--mu-r`` so the panel's preset label
        (e.g. "Steel (mu_r=100)") doesn't need to map 1:1 to the
        add_material_args choices (steel/copper/aluminum/custom).

        CLI-DIFF: ignore --reg --shift-eps --nthreads --output -- advanced solver knobs and auto-generated output paths; deliberately defaulted.
        CLI-DIFF: ignore --output -- auto-injected by calc_main wrapper (calc_common.py).
        """
        step = self.val("peec_step")
        if not step:
            raise ValueError("PEEC STEP file is required for "
                              "PEEC+FEM+Kelvin.")
        if not os.path.isfile(step):
            raise ValueError(f"STEP file not found: {step}")
        solver = self._FEM_SOLVER_MAP.get(self.val("solver"), "pardiso")
        # The method name promises "+ Kelvin" so enforce that the .vol
        # actually has a kelvin material with periodic ID -- otherwise
        # calc_fem_kelvin silently downgrades to reg-only gauge
        # truncation (RISK 3 fix, 2026-05-12).
        cmd = [_PYTHON, calc_script("calc_fem_kelvin.py"),
               "--vol", vol_path,
               "--fes-order", str(self.val("fes_order")),
               "--frequency", self.val("freq"),
               "--material", "custom",
               "--sigma", self.val("wp_sigma"),
               "--mu-r", self.val("mu_r"),
               "--impedance", self._impedance_model_cli(),
               "--formulation", "total",
               "--current", self.val("current"),
               "--half-thickness", self.val("half_thickness"),
               "--solver", solver,
               "--peec-step", step,
               "--peec-sigma", self.val("coil_sigma"),
               "--peec-nwinc", str(self.val("peec_nwinc")),
               "--peec-nhinc", str(self.val("peec_nhinc")),
               "--peec-n-peri", str(self.val("peec_n_peri")),
               "--require-kelvin",
               "--msh-output", msh_output(vol_path, "_fem_kelvin"),
               "--output", json_output(vol_path, "_fem_kelvin")]
        if self._impedance_model_cli() == "esim":
            bh = self.val("bh_file")
            if not bh:
                raise ValueError(
                    "Nonlinear ESIM impedance model requires a BH "
                    "file (2-column H[A/m] B[T]).  Fill the 'BH file' "
                    "field, or switch to Linear SIBC.")
            cmd += ["--bh-file", bh,
                    "--max-iter", str(self.val("esim_max_iter"))]
            # calc_fem_kelvin uses --max-iter, not --esim-max-iter.
            # --esim-tol / --esim-relax / --esim-anderson-m have no
            # equivalent (verified 2026-05-24 vs calc_fem_kelvin.py
            # argparse).  Per-DOF ESIM IS supported via
            # --esim-per-panel.
            if self._widgets["esim_per_panel"].isChecked():
                cmd.append("--esim-per-panel")
        return cmd

    def _build_fem_coilmesh_command(self, vol_path):
        solver = self._FEM_SOLVER_MAP.get(self.val("solver"), "pardiso")
        # METHOD_FEM_FULL's label is "Full simulation (FEM A-V + wp SIBC
        # + Kelvin)" -- enforce that the .vol actually carries the
        # Kelvin extension (RISK 3 fix, 2026-05-12).
        cmd = [_PYTHON, calc_script("calc_fem_coilmesh.py"),
               "--vol", vol_path,
               "--frequency", self.val("freq"),
               "--current", self.val("current"),
               "--coil-sigma", self.val("coil_sigma"),
               "--sigma", self.val("wp_sigma"),
               "--mu-r", self.val("mu_r"),
               "--half-thickness", self.val("half_thickness"),
               "--fes-order", str(self.val("fes_order")),
               "--solver", solver,
               "--sibc-bnd", "sibc",
               "--source-bnd", "source",
               "--sink-bnd", "sink",
               "--coil-mat", "coil",
               "--impedance-model", self._impedance_model_cli(),
               "--require-kelvin",
               "--msh-output", msh_output(vol_path, "_fem_full"),
               "--output", json_output(vol_path, "_fem_full")]
        if self._impedance_model_cli() == "esim":
            bh = self.val("bh_file")
            if not bh:
                raise ValueError(
                    "Nonlinear ESIM impedance model requires a BH "
                    "file (2-column H[A/m] B[T]).  Fill the 'BH file' "
                    "field, or switch to Linear SIBC.")
            cmd += ["--bh-file", bh,
                    "--esim-max-iter", str(self.val("esim_max_iter")),
                    "--esim-tol", self.val("esim_tol"),
                    "--esim-relax", self.val("esim_relax")]
            # calc_fem_coilmesh.py accepts --esim-per-panel + --esim-relax
            # but NOT --esim-anderson-m (verified 2026-05-24 vs argparse).
            if self._widgets["esim_per_panel"].isChecked():
                cmd.append("--esim-per-panel")
        return cmd


# ============================================================
# IH Window (with .vol-load-hook + formatted output)
# ============================================================

class IHWindow(AnalysisWindow):
    def __init__(self, vol_path=""):
        super().__init__("Radia - Induction Heating", vol_path,
                         settings_key="ih")
        panel = IHPanel()
        self._set_panel(panel)
        # Restore previous session's widget state FIRST, then let any
        # explicit constructor-supplied vol_path override the restored
        # wp_vol.  The previous order (setText before restore) silently
        # lost the Cubit Solve menu's current .vol path whenever a
        # different .vol was saved last session, because restore would
        # overwrite the just-set wp_vol text.
        self._restore_settings()
        if vol_path and "wp_vol" in panel._widgets:
            panel._widgets["wp_vol"].setText(self.display_path(vol_path))
        # Re-fire method change so visibility hooks run after the panel-
        # window binding (wp_vol section visibility) is wired.
        panel._on_method_changed(panel.val("method"))
        # Trigger label re-inspection for the panel-owned wp_vol widget
        # (in case the constructor / settings restored a non-empty path).
        if "wp_vol" in panel._widgets:
            panel._on_wp_vol_changed_text(panel._widgets["wp_vol"].text())
        # PEEC-inductance convenience: if the method is PEEC-inductance
        # AND the STEP field is empty (first launch / fresh QSettings),
        # auto-populate from the newest *.step / *.stp in cwd.
        # This was previously a separate PEECInductanceWindow class;
        # merged into IHWindow 2026-04-26 because PEEC-inductance is
        # just one of IH's four methods and a standalone window only
        # added a thin wrapper that forced the method.
        self._maybe_auto_fill_step_from_cwd()
        # "Run thermal..." chain button.  Inserted next to Open GMSH
        # in the action row.  Enabled when an IH solve emits a usable
        # qsurf.sol companion (Phase A keys "qsurf_sol" / "msh_file");
        # click switches the Method dropdown to Thermal (the embedded
        # HeatPanel from _heat_panel.py) with the EM-side outputs
        # pre-filled so the user picks only the workpiece thermal
        # mesh and the thermal parameters.
        self._heat_qsurf_sol = ""
        self._heat_em_vol = ""
        self._heat_btn = self._install_heat_button()
        self._update_run_state()

    def _install_heat_button(self):
        """Append a 'Run thermal...' button to the AnalysisWindow
        action row (Run / Stop / Open GMSH).  Initially disabled."""
        from PySide6.QtWidgets import QPushButton, QStyle
        btn_row = self._gmsh_btn.parent().layout()
        if btn_row is None:
            return None
        style = self.style()
        btn = QPushButton(
            style.standardIcon(QStyle.SP_ArrowRight), " Switch to Thermal")
        btn.setFixedHeight(32)
        btn.setEnabled(False)
        btn.setToolTip(
            "Switch the Method dropdown to 'Thermal: 3D + rotation' so "
            "you can configure a heat-transfer run.<br><br>"
            "<b>Does NOT auto-fill</b> any fields -- the .sol + .vol "
            "from the just-completed EM run are PRINTED to the output "
            "box so you can copy-paste them into the Thermal section's "
            "<i>qsurf .sol</i> and <i>EM .vol</i> Browse pickers (or "
            "pick different files entirely).  Explicit-input is by "
            "design (2026-05-24) -- pre-2026-05-24 magic auto-fill "
            "led to silent wrong-pair selection.")
        btn.clicked.connect(self._on_run_thermal)
        gmsh_idx = btn_row.indexOf(self._gmsh_btn)
        btn_row.insertWidget(gmsh_idx + 1, btn)
        return btn

    def _maybe_auto_fill_step_from_cwd(self):
        """If the active method needs a STEP coil AND the field is
        empty, populate it with the newest .step / .stp in the current
        working directory.  No-op when the user already has a value.
        Layer 3 panels never accept .jou (Cubit-only Layer 1 input);
        the calc layer is fed STEP / .vol only."""
        panel = self._panel
        if panel is None:
            return
        method = panel.val("method")
        # Only PEEC-inductance auto-fills; other PEEC methods (BEM /
        # FEM-Kelvin) ALSO use peec_step but the wp .vol is the user's
        # primary intent so we do not auto-fill them.
        if method != METHOD_PEEC_IND:
            return
        step_widget = panel._widgets.get("peec_step")
        if step_widget is None or step_widget.text().strip():
            return
        import glob
        candidates = sorted(
            glob.glob("*.step") + glob.glob("*.stp"),
            key=lambda p: os.path.getmtime(p), reverse=True)
        if candidates:
            step_widget.setText(os.path.abspath(candidates[0]))

    def _on_finished(self, exit_code, exit_status):
        # Delegate core finish handling + then append IH-specific summary.
        super()._on_finished(exit_code, exit_status)
        # Reset chain state on every run -- "Run thermal..." is only
        # offered for the most recent successful EM solve.
        self._heat_qsurf_sol = ""
        self._heat_em_vol = ""
        if self._heat_btn is not None:
            self._heat_btn.setEnabled(False)
        if exit_code != 0:
            return
        # Parse JSON one more time for IH-specific pretty-print
        try:
            import json
            text = self._output.toPlainText()
            result = None
            for line in reversed(text.split("\n")):
                line = line.strip()
                if line.startswith("{"):
                    try:
                        result = json.loads(line)
                        break
                    except json.JSONDecodeError:
                        pass
            if result is None or "error" in result:
                return

            # Detect a usable qsurf.sol pair from the JSON result.
            # Phase A's calc_fem_kelvin emits "qsurf_sol" directly;
            # the companion EM .vol is at "<msh stem>_fem.vol" by
            # the save_vol_sol_pair convention.  Other IH methods
            # (PEEC inductance / PEEC+BEM) do not save qsurf.sol so
            # the chain button stays disabled for them.
            qsurf = result.get("qsurf_sol") or ""
            if qsurf and os.path.isfile(qsurf):
                em_vol = ""
                msh_file = result.get("msh_file") or ""
                if msh_file:
                    msh_stem = os.path.splitext(msh_file)[0]
                    candidate = msh_stem + "_fem.vol"
                    if os.path.isfile(candidate):
                        em_vol = candidate
                if not em_vol and qsurf.endswith("_qsurf.sol"):
                    candidate = qsurf[:-len("_qsurf.sol")] + "_fem.vol"
                    if os.path.isfile(candidate):
                        em_vol = candidate
                # Final fallback: try the wp_vol's directory + the JSON's
                # ``wp_vol`` (an additional convention some IH methods
                # emit when the .vol carries the EM mesh directly).
                if not em_vol:
                    wp_v = result.get("wp_vol") or result.get("vol") or ""
                    if wp_v and os.path.isfile(wp_v):
                        em_vol = wp_v
                self._heat_qsurf_sol = qsurf
                self._heat_em_vol = em_vol
                if self._heat_btn is not None:
                    self._heat_btn.setEnabled(True)
                if not em_vol:
                    # v4.58.0+ contract: thermal panel's Run button stays
                    # disabled until the user fills --em-vol explicitly
                    # (.sol is a coefficient vector only, NGSolve cannot
                    # load it without the matching mesh).  Surface a
                    # hint so the user knows to pick the .vol manually
                    # in the thermal window after launch.
                    self._output.appendPlainText(
                        "\n[hint] EM .vol companion to qsurf.sol was not "
                        "auto-located.  When the thermal panel opens, "
                        "you'll need to browse to the EM .vol "
                        "(usually <stem>_fem.vol next to <stem>_qsurf.sol) "
                        "before the Run button enables.")

            lines = ["", "=== IH Summary ==="]
            method = result.get("method", "")
            if method:
                lines.append(f"  Method: {method}")
            # Inductance
            for key, label in (("L_coil_nH", "L_coil (vacuum)"),
                                ("L_total_nH", "L_total (with wp)")):
                if key in result:
                    lines.append(f"  {label}: {result[key]:.3f} nH")
            if "delta_L_nH" in result and result["delta_L_nH"] is not None:
                tag = ""
                rel = result.get("delta_L_reliability", "")
                if rel == "experimental":
                    tag = " [EXPERIMENTAL]"
                form_label = "Telegen φ·B" if rel == "production" else "Telegen"
                lines.append(
                    f"  ΔL_wp ({form_label}): {result['delta_L_nH']:+.3f} nH"
                    f"{tag}")
                if "delta_L_JsA_nH" in result \
                        and result["delta_L_JsA_nH"] is not None:
                    lines.append(
                        f"    diag J_s·A:  "
                        f"{result['delta_L_JsA_nH']:+.3f} nH "
                        f"(continuum-equivalent)")
            if "R_coil_mOhm" in result:
                lines.append(
                    f"  R_coil: {result['R_coil_mOhm']:.4f} mOhm")
            if "R_total_mOhm" in result and result["R_total_mOhm"] is not None:
                lines.append(
                    f"  R_total: {result['R_total_mOhm']:.4f} mOhm "
                    f"(coil + wp Telegen)")
            if ("delta_R_mOhm" in result and result["delta_R_mOhm"] is not None
                and "delta_R_expected_mOhm" in result
                and result["delta_R_expected_mOhm"] is not None):
                dR = result["delta_R_mOhm"]
                dRe = result["delta_R_expected_mOhm"]
                ratio = (dR / dRe) if abs(dRe) > 1e-30 else 0.0
                tag = "" if 0.5 < ratio < 2.0 else \
                    f" [WARN: {ratio:.3f}x energy-balance]"
                lines.append(
                    f"  ΔR_wp (Telegen): {dR:+.4f} mOhm  "
                    f"(energy-balance: {dRe:+.4f}){tag}")
            if "n_filaments" in result:
                lines.append(f"  filaments: {result['n_filaments']}")
            # Dissipation
            if "P_wp_W" in result:
                p = result["P_wp_W"]
                lines.append(f"  P_workpiece: {p:.4e} W  ({p*1e3:.3f} mW)")
            if "P_coil_W" in result:
                p = result["P_coil_W"]
                lines.append(f"  P_coil (ohmic):    {p:.4e} W "
                              f"(mesh-sensitive ±15%)")
            if "P_total_W" in result:
                p = result["P_total_W"]
                lines.append(f"  P_total: {p:.4e} W")
            if "P_wp_W" in result and "P_total_W" in result and \
                    result["P_total_W"] > 0:
                eta = result["P_wp_W"] / result["P_total_W"] * 100
                lines.append(f"  Heating efficiency: {eta:.1f}% "
                              f"(P_wp/P_total)")
            # H_t / area
            if "H_t_rms_wp_Am" in result:
                lines.append(f"  H_t_rms (wp): "
                              f"{result['H_t_rms_wp_Am']:.2f} A/m")
            if "wp_area_m2" in result:
                lines.append(
                    f"  wp area: {result['wp_area_m2']*1e4:.2f} cm^2")
            # Coil diagnostics
            if "coil_delta_mm" in result:
                lines.append(
                    f"  coil delta = {result['coil_delta_mm']:.3f} mm")
            if "coil_h_max_mm" in result and "coil_delta_mm" in result:
                h = result["coil_h_max_mm"]
                d = result["coil_delta_mm"]
                tag = " [OK]" if h <= d else " [WARN under-resolved]"
                lines.append(
                    f"  coil h_max = {h:.3f} mm "
                    f"(delta = {d:.3f} mm){tag}")
            # Timings
            if "t_solve_s" in result:
                lines.append(f"  solve: {result['t_solve_s']:.1f} s")

            # ------ Thermal-method summary ------
            # calc_heat / calc_heat_axisym emit T_max_C / T_min_C /
            # Q_input_J / probe history / msh_file / vtu_files.  No
            # L/R/P keys (those are EM-side).  Show the temperature
            # stats so the user does not need to open the .msh just
            # to see whether the workpiece reached soak temperature.
            t_max = result.get("T_max_C")
            t_min = result.get("T_min_C")
            if t_max is not None or t_min is not None:
                lines.append("  --- Thermal ---")
            if t_max is not None:
                lines.append(f"  T_max: {t_max:.2f} degC")
            if t_min is not None:
                lines.append(f"  T_min: {t_min:.2f} degC")
            t_init = result.get("T_initial_C")
            if t_init is not None and t_max is not None:
                rise = t_max - t_init
                lines.append(f"  Delta T (T_max - T_initial): {rise:+.2f} K")
            qin = result.get("Q_input_J")
            if qin is not None:
                lines.append(f"  Q_input: {qin:.3e} J  ({qin/1e3:.3f} kJ)")
            t_end = result.get("t_end_s")
            n_steps = result.get("n_steps")
            if t_end is not None and n_steps is not None:
                lines.append(
                    f"  Time: {t_end:.2f} s in {int(n_steps)} steps")
            # Rotation status (v4.58.0+).
            rpm = result.get("rotation_rpm")
            if rpm is not None and float(rpm) > 0:
                lines.append(
                    f"  Rotation: {float(rpm):g} rpm "
                    f"(spinning workpiece)")
            # Probe at the user's chosen point (if any) -- show final
            # value and rise.  The full history is in T_probe_history_C
            # which is too long for inline display.
            probe_hist = result.get("T_probe_history_C")
            if probe_hist and isinstance(probe_hist, list) \
                    and len(probe_hist) >= 2:
                t_probe_final = probe_hist[-1]
                t_probe_initial = probe_hist[0]
                lines.append(
                    f"  Probe: {t_probe_initial:.2f} -> "
                    f"{t_probe_final:.2f} degC "
                    f"({t_probe_final - t_probe_initial:+.2f} K)")
            # Output file paths -- so the user can see at a glance
            # what got written.
            t_sol = result.get("T_sol_file") or ""
            heat_vol = result.get("heat_vol_file") or ""
            if t_sol:
                lines.append(
                    f"  T .sol: {os.path.basename(t_sol)}  "
                    f"(re-loadable for post-processing)")
            if heat_vol:
                lines.append(
                    f"  heat .vol: {os.path.basename(heat_vol)}  "
                    f"(companion mesh for T .sol)")
            msh = result.get("msh_file") or ""
            if msh:
                lines.append(f"  GMSH .msh: {os.path.basename(msh)}")
            vtu_files = result.get("vtu_files") or []
            if vtu_files:
                lines.append(
                    f"  VTU files: {len(vtu_files)} steps "
                    f"({os.path.basename(vtu_files[0])} ... "
                    f"{os.path.basename(vtu_files[-1])})")

            self._output.appendPlainText("\n".join(lines))

            # Auto-fire "Open GMSH" for the Thermal method only.  EM
            # methods leave the button enabled but do not auto-open,
            # matching the prior UX (Thermal is the new "closing
            # step" of the IH pipeline so the T distribution should
            # appear in front of the user without an extra click).
            if (result.get("method") == "thermal-3d"
                    or result.get("method") == "thermal-axisym"
                    or (t_max is not None and msh)):
                if self._last_msh and self._gmsh_btn \
                        and self._gmsh_btn.isEnabled():
                    self._output.appendPlainText(
                        "\n(Auto-opening GMSH on the T distribution...)")
                    self._open_gmsh()
        except Exception as e:
            self._output.appendPlainText(f"(IH summary skipped: {e})")

        # Re-persist the .log AFTER the IH-specific summary block.  The
        # base class's super()._on_finished call wrote a .log that
        # ended at "--- Result ---" + the generic _append_standard_summary
        # block.  All the IH-specific lines above (L_coil, P_workpiece,
        # T_max, file paths, auto-open banner) landed in _output AFTER
        # that, so they were missing from the .log -- which defeats the
        # Result Output Persistence Policy's triage purpose.  Overwrite
        # the .log now so it captures the complete on-screen output.
        self._persist_output_log()

    def _on_run_thermal(self):
        """Switch the method dropdown to Thermal -- explicit-input UX.

        Pre-2026-05-24 this button auto-filled ``qsurf_sol`` and
        ``em_vol`` from the just-completed EM solve.  The auto-fill
        was confusing in practice: users landed on the Thermal panel
        with pre-populated fields and could not tell whether the
        values were from this session's EM run, a QSettings restore,
        or something else.  Two adjacent EM runs that wrote different
        ``.sol`` files would silently inherit the wrong pair if the
        user didn't notice the field contents.

        Per user feedback 2026-05-24 ("Auto-Fillはわかりにくいので陽に
        指定する形に変更しよう"), the button now does ONLY the method
        switch.  Browse-selecting the .sol + .vol is the user's
        responsibility and the file paths are visible in the Browse
        widgets, so no surprise inheritance is possible.

        Heat analysis still lives in the same IHWindow under the
        Thermal method choice; this button is just a shortcut for
        the method dropdown change (one click vs digging through the
        list).
        """
        panel = self._panel
        heat = getattr(panel, "_heat_panel", None)
        if heat is None:
            self._output.appendPlainText(
                "\nThermal sub-panel not available (_heat_panel.HeatPanel "
                "import failed at startup).")
            return
        # Switch the method dropdown to the most general Thermal
        # entry (3D + rotation -- handles both rotating and static
        # cases since rotation_rpm=0 yields the static behaviour).
        # The user can pick a different Thermal method (axisym /
        # 3D static) afterwards if they prefer.
        from radia.radia_ih import METHOD_THERMAL_3D_ROTATING
        panel._method_combo.setCurrentText(METHOD_THERMAL_3D_ROTATING)
        self._output.appendPlainText(
            f"\nSwitched to '{METHOD_THERMAL_3D_ROTATING}'.")
        # Surface the .sol / .vol paths from the EM run as TEXT ONLY
        # (not auto-populated into the Browse fields).  User copies
        # them into the Browse pickers if they want this pair, OR
        # picks a different qsurf .sol from elsewhere.  This keeps
        # provenance visible without silently inheriting state.
        if self._heat_qsurf_sol:
            self._output.appendPlainText(
                f"  Last EM run wrote: {self._heat_qsurf_sol}")
            if self._heat_em_vol:
                self._output.appendPlainText(
                    f"                  +  {self._heat_em_vol}")
            self._output.appendPlainText(
                "  Browse qsurf_sol + em_vol in the Thermal section "
                "to use these (the panel does NOT auto-fill them).")
        else:
            self._output.appendPlainText(
                "  No qsurf .sol from this session -- Browse manually "
                "in the Thermal section.")


def main():
    run_app(IHWindow)


if __name__ == "__main__":
    main()
