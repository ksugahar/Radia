"""Radia IH (Induction Heating) analysis window.

Two methods (2026-04-19):

  Fast workpiece heating (PEEC+BEM, 1-way)
    -> calc_peec_bem.py   (~3 min, P_wp ±5%, L_coil vacuum only)

  Full simulation (FEM A-V + wp SIBC + Kelvin)
    -> calc_fem_coilmesh.py (~1-7 min, L + P_wp + P_coil)

Both drive from a gapped torus coil (real IH has physical port
terminations — closed-torus topology is not supported).
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
             METHOD_PEEC_FEM_KELVIN, METHOD_FEM_FULL])
        self._method_combo.currentTextChanged.connect(self._on_method_changed)
        self._method_combo.setToolTip(METHOD_TOOLTIP[METHOD_PEEC_IND])

        # Per-method status line (.vol label check / skin depth hint).
        # Base helper: idiom for de-emphasised status row; updates via
        # self._status_label.setText(...) elsewhere in this class.
        self.add_status_label()

        # ============ Drive (frequency + current) ============
        self._add_section("Drive")
        freq = self.add_line("freq", "Frequency [Hz]:", "7000")
        freq.editingFinished.connect(self._update_status)
        freq.setToolTip("Operating frequency.  Skin depth\n"
                         "  delta = sqrt(2 / (w mu sigma))")

        current = self.add_line("current", "Coil current [A, peak]:", "1.0")
        current.setToolTip(
            "Peak (not RMS) coil port current.  Field quantities are "
            "complex phasors; output P is time-averaged (1/2 Re).")

        # ============ Coil material (INDEPENDENT from WP) ============
        self._add_section("Coil material")
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
             "Nonlinear ESIM (experimental, WIP)"],
            default=0)
        imp.currentTextChanged.connect(self._on_impedance_changed)
        imp.setToolTip(
            "<b>Linear SIBC</b>: Z_s = (1+j) rho/delta * sqrt(mu_r). "
            "Ok for Cu/Al, and for steel with a constant mu_r.<br>"
            "<b>ESIM</b>: 1D cell problem solves B-H(H) self-consistently "
            "(Karl iteration). Needed when mu_r varies with H (saturated "
            "steel). <b>Calc script support is WIP</b> — panel accepts "
            "settings but the subprocess may reject.")

        # ESIM-only widgets
        self.add_browse("bh_file", "BH file:", default="",
                         filter_str="BH tables (*.txt *.csv);;All (*)")
        self.add_spin("esim_max_iter", "max iter:", 15, 1, 200)
        self.add_line("esim_tol", "tolerance:", "1e-3")

        # ============ Linear solver (method-dependent) ============
        self._add_section("Linear solver")
        solver = self.add_combo("solver", "Solver:", ["pardiso"])
        solver.setToolTip(
            "<b>PEEC+BEM</b>:<br>"
            "  Dense LU — fast for small filament bundles (<500 segments)<br>"
            "  HACApK — O(N log N), for large bundles<br>"
            "<br>"
            "<b>FEM A-V</b>:<br>"
            "  pardiso — sparse direct (default, fast, memory-heavy)<br>"
            "  AMS — Compact AMS+COCR for HCurl p=1 (low memory; shifted preconditioner internally)<br>"
            "  BDDC — preconditioned CG, recommended for p&gt;=2<br>"
            "  iccg — generic fallback (Incomplete Cholesky + CG)")

        # ============ Advanced (collapsed by default) ============
        self._add_section("Advanced")
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
        # time via ``radia_export netgen "f.vol" order N``) and
        # auto-detected by calc_inductance.py from the companion
        # ``.vol.json``.  A post-load ``mesh.Curve(p)`` silently falls
        # back to flat without a CAD callback, so we never expose it.
        self.add_spin("fes_order", "Basis order:", 1, 1, 3)

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
        for key in ("bh_file", "esim_max_iter", "esim_tol"):
            self._set_row_visible(key, is_esim)
        self._update_status()

    # ----------------------- Method + visibility -----------------------

    def _on_method_changed(self, method):
        is_peec_ind = (method == METHOD_PEEC_IND)
        is_bema_ind = (method == METHOD_BEMA_IND)
        is_peec_bem = (method == METHOD_PEEC_BEM)
        is_bema_bem = (method == METHOD_BEMA_BEM)
        is_peec_fem_k = (method == METHOD_PEEC_FEM_KELVIN)
        is_fem = (method == METHOD_FEM_FULL)
        # Inductance-style (vacuum) vs weak-coupled (workpiece) groupings.
        is_vacuum = is_peec_ind or is_bema_ind
        is_weak = is_peec_bem or is_bema_bem
        is_inductance_path = is_vacuum or is_weak    # any calc_inductance.py mode
        is_bem_a = is_bema_ind or is_bema_bem
        is_peec = is_peec_ind or is_peec_bem
        # PEEC + PEEC+FEM+Kelvin take CAD STEP (filament coil); BEM-A
        # variants take a pre-meshed .vol coil.  FEM-full uses a
        # volumetric coil baked into its workpiece .vol.
        needs_step = is_peec or is_peec_fem_k
        needs_coil_vol = is_bem_a
        needs_wp = is_weak or is_peec_fem_k or is_fem

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
        fes_spin.setMaximum(2 if is_weak else 3)
        # Workpiece geometry (.vol) is shown for all modes that need a
        # workpiece mesh; the legacy top-level Model row is gone.
        self._set_row_visible("_sec_wp_vol", needs_wp)
        self._set_row_visible("wp_vol", needs_wp)
        self._set_row_visible("_sec_wp_material", needs_wp)
        self._set_row_visible("_sec_wp_imp", needs_wp)

        # Workpiece widgets are meaningless for vacuum-only modes.
        for key in ("wp_material", "wp_sigma", "mu_r", "half_thickness",
                    "impedance_model", "bh_file",
                    "esim_max_iter", "esim_tol"):
            self._set_row_visible(key, needs_wp)

        # ESIM sub-widgets re-evaluated by _on_impedance_changed when
        # needs_wp is True.
        if needs_wp:
            self._on_impedance_changed(self.val("impedance_model"))

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

        self._status_label.setText("<br>".join(lines))

    # ----------------------- External hooks -----------------------

    def is_runnable(self):
        method = self._method_combo.currentText()
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
        if "wp_vol" not in self._widgets:
            return ""
        return self.val("wp_vol")

    def _on_wp_vol_changed_text(self, _text):
        """Re-inspect the workpiece .vol labels after the user edits the
        wp_vol field (Browse... or manual edit).  Mirrors the legacy
        IHWindow._on_vol_changed wiring against the panel-owned widget.
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
            coil_arg = ["--coil-step", coil_in]
        else:  # bem-a
            coil_in = self.val("coil_vol") if "coil_vol" in self._widgets else ""
            if not coil_in:
                raise ValueError("Coil .vol file is required for BEM-A weak-coupled mode.")
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
            if bh:
                cmd += ["--bh-file", bh,
                        "--esim-max-iter", str(self.val("esim_max_iter")),
                        "--esim-tol", self.val("esim_tol")]
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
               "--msh-output", msh_output(vol_path, "_fem_kelvin"),
               "--output", json_output(vol_path, "_fem_kelvin")]
        if self._impedance_model_cli() == "esim":
            bh = self.val("bh_file")
            if bh:
                cmd += ["--bh-file", bh,
                        "--max-iter", str(self.val("esim_max_iter"))]
                # calc_fem_kelvin uses --max-iter, not --esim-max-iter
                # --esim-tol has no equivalent; ESIM tolerance is hardcoded
        return cmd

    def _build_fem_coilmesh_command(self, vol_path):
        solver = self._FEM_SOLVER_MAP.get(self.val("solver"), "pardiso")
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
               "--msh-output", msh_output(vol_path, "_fem_full"),
               "--output", json_output(vol_path, "_fem_full")]
        if self._impedance_model_cli() == "esim":
            bh = self.val("bh_file")
            if bh:
                cmd += ["--bh-file", bh,
                        "--esim-max-iter", str(self.val("esim_max_iter")),
                        "--esim-tol", self.val("esim_tol")]
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
        # click launches radia_heat.py with the EM-side outputs
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
            style.standardIcon(QStyle.SP_ArrowRight), " Run thermal...")
        btn.setFixedHeight(32)
        btn.setEnabled(False)
        btn.setToolTip(
            "Launch the thermal panel (radia_heat) with this run's "
            "q_surf .sol pre-filled.  Pick a workpiece-volume .vol "
            "in the new window; the EM-side q_surf is projected "
            "onto its heating surface automatically.")
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
                self._heat_qsurf_sol = qsurf
                self._heat_em_vol = em_vol
                if self._heat_btn is not None:
                    self._heat_btn.setEnabled(True)

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

            self._output.appendPlainText("\n".join(lines))
        except Exception as e:
            self._output.appendPlainText(f"(IH summary skipped: {e})")

    def _on_run_thermal(self):
        """Launch radia_heat.py with this run's qsurf .sol pre-filled.

        We start a detached subprocess (CREATE_NEW_PROCESS_GROUP on
        Windows so closing the IH window doesn't kill the thermal
        window).  No state is shared with the IH window after launch
        beyond the CLI arguments below.
        """
        if not self._heat_qsurf_sol:
            return
        import subprocess
        heat_script = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "radia_heat.py")
        cmd = [_PYTHON, heat_script,
               "--qsurf-sol", self._heat_qsurf_sol]
        if self._heat_em_vol:
            cmd += ["--em-vol", self._heat_em_vol]
        # CREATE_NEW_PROCESS_GROUP = 0x00000200 lets the child
        # outlive the IH parent without needing fork().
        flags = 0x00000200 if sys.platform == "win32" else 0
        try:
            subprocess.Popen(cmd, creationflags=flags)
            self._output.appendPlainText(
                f"\nLaunched thermal panel: {os.path.basename(heat_script)}")
            self._output.appendPlainText(
                f"  qsurf_sol = {self._heat_qsurf_sol}")
            if self._heat_em_vol:
                self._output.appendPlainText(
                    f"  em_vol    = {self._heat_em_vol}")
        except Exception as e:
            self._output.appendPlainText(
                f"\nFailed to launch thermal panel: {e}")


def main():
    run_app(IHWindow)


if __name__ == "__main__":
    main()
