"""
Radia GUI / Cubit panel pitfalls — bugs we hit, root causes, and the rules
that prevent them from coming back. Read this BEFORE adding new
parameters, methods, or modes to a `radia_*.py` panel.

Each entry is structured as:

  ## <symptom>
  Root cause: ...
  Rule: ...
  Reference: <commit hash, file:line>

The MCP server exposes this via panel_gui_pitfalls(topic=...). Topics
are stable keywords (combo_state, mode_switch, gmsh_viz, subprocess_args,
cubit_jou, wheel_guard, panel_layout, result_output) so they can be
referenced from prompts.
"""

PANEL_GUI_PITFALLS = """\
# Radia GUI / Cubit Panel Pitfalls

This document lists pitfalls we have hit (and fixed) while developing
the Radia panels and the calc_*.py subprocess scripts. Each section
ends with a **rule** to prevent the same bug from coming back.

Topics: combo_state, mode_switch, gmsh_viz, subprocess_args, cubit_jou,
        sample_jou, layout_unification, learn_edition_cap,
        wheel_guard, panel_layout, result_output

============================================================
## combo_state — Combo box save/restore by INDEX is fragile
============================================================

**Symptom**: After dropping a Method combo item ("BEM-SIBC (WP)") in
one release, the next time the panel opens the Method combo is BLANK.
The user cannot pick anything until they manually re-select.

**Root cause**: ``ModePanel.save_state`` was storing
``QComboBox.currentIndex()`` (an integer). The user's saved settings
contained ``method=2`` from the previous 3-item combo
(BEM, BEM-SIBC (WP), FEM). After the cleanup the combo had only 2
items (BEM, FEM). ``setCurrentIndex(2)`` on a 2-item combo silently
sets the selection to ``-1`` (no item displayed).

**Rule**: Save QComboBox state by **text** (``currentText()``), not
index. On restore use ``findText(val)`` and only apply when
``>= 0``; otherwise leave the panel default. Index-based restore is
acceptable ONLY as a legacy migration path with bounds-check
(``0 <= ival < w.count()``).

**Reference**: commit ``1aa66ee``, ``radia_gui_base.py::ModePanel.
save_state / restore_state``.

============================================================
## mode_switch — Hidden widgets still feed build_command
============================================================

**Symptom**: User picks FEM mode, runs solve, the script crashes with
``argument --material: invalid choice: 'custom'`` even though the FEM
panel does NOT show any custom-material UI. Or: the FEM solve uses
``mu_r = 100`` (the BEM default) regardless of what the user picks.

**Root cause**: ``_build_fem_command`` reads ``self.val("mu_r")``,
but the ``mu_r`` widget was created in the BEM workpiece group and
hidden in FEM mode. The user has no way to change it from the FEM
panel, so the BEM default leaks into the FEM call. Worse, the
``--material custom`` flag was sent without the calc script declaring
"custom" as an accepted ``argparse`` choice — the script died with
exit 2.

**Rules**:

  1. Any widget READ by ``_build_*_command`` for a method MUST be
     visible in that method's panel. Hidden widgets are
     no-user-control state — silent UI bug.
  2. Every CLI flag the GUI sends MUST be in the calc script's
     ``argparse choices``. Use ``custom`` as the explicit "user
     specifies mu_r/sigma directly" choice when the script also
     supports built-in material lookup.
  3. When you delete a Method combo entry, scrub ALL widget
     visibility in ``_on_method_changed`` so the new mode does not
     inherit stale rows from the old one.

**Reference**: commit ``8b34b04`` (calc_fem_kelvin --material custom),
commit ``1aa66ee`` (unified BEM/FEM widget set).

============================================================
## layout_unification — Mode-specific widgets are usually wrong
============================================================

**Symptom**: BEM panel shows "Coil sigma" but FEM does not. FEM shows
"Current [A]" + "Coil radius" + "Workpiece radius" but BEM does not.
The two panels look unrelated even though they solve the same
physics problem.

**Root cause**: Each method was developed independently with its own
widget set. Over time the two diverged because nobody enforced
"same physics -> same parameters".

**Rule**: For panels that offer multiple solver methods, the widget
set should be a SHARED superset. Method-specific differences should
be limited to:

  - Solver combo items (different solver per method)
  - Iteration cap (FEM iterative, BEM direct)
  - Post-processing toggles (e.g. "air field calc" — BEM only because
    FEM volume mesh always carries the field)

For Radia IH the user has stated explicitly:

  > BEMとFEMの違いは、Solverの違いと、Air field calcの有無程度では？

i.e. the only allowed differences are Solver and Air field calc.
Frequency, current, coil sigma, workpiece sigma, mu_r, half thickness,
SIBC/ESIM choice, BH file are SHARED across methods.

**Reference**: commit ``1aa66ee`` from the retired desktop-panel era
("Layout (per user request 2026-04-12)" comment block).  Current panels
should express the shared settings through `IHDesignSpec` and the notebook
workbench, not by recreating desktop widgets.

============================================================
## gmsh_viz — Mesh elements obscure the field
============================================================

**Symptom**: GMSH opens the .msh file the panel wrote and the model
is **completely black** — a sphere of internal triangle edges
covering the entire object. No field arrows visible.

**Root cause**: ``GmshPostExport`` writes the volume / surface mesh
together with the field views. GMSH defaults render
``Mesh.Volumes / VolumeFaces / VolumeEdges / SurfaceEdges /
SurfaceFaces`` all to ON, so a 100k-element sphere mesh draws as a
black blob.

**Rule**: Always write a **companion ``.geo`` file** next to the
``.msh`` that hides the mesh primitives and configures the field
views as 3D arrows:

```
Merge "model.msh";
Mesh.NumSubEdges = 4;
Mesh.Volumes = 0;
Mesh.VolumeFaces = 0;
Mesh.VolumeEdges = 0;
Mesh.SurfaceEdges = 0;
Mesh.SurfaceFaces = 0;
View[0].VectorType = 4;       // 3D arrow
View[0].IntervalsType = 3;
View[1].VectorType = 4;
View[1].IntervalsType = 3;
```

The "Open GMSH" button in the panel should open the ``.geo``, NOT
the raw ``.msh``, so these display options are applied automatically.

**Reference**: commit ``1aa66ee``, ``calc_fem_kelvin.py`` Step 7;
commit ``8dbee35``, ``calc_inductance.py::post_process``.

============================================================
## gmsh_viz — Scalar |B| / |J| views are redundant
============================================================

**Symptom**: The panel exports both a vector ``B`` view and a scalar
``|B|`` view. The user reports they only want the vector. Same for
``|J|``.

**Root cause**: GMSH's vector view already shows the magnitude (the
arrow color encodes ``|B|`` by default). A separate scalar view
clutters the View list and wastes file space.

**Rule**: Export VECTORS only. Drop ``add_field("|B|", ...)`` and
``add_scalar_field("|B|", ...)`` calls. The user picks the vector
view in GMSH and gets magnitude as a free side effect. The same
applies to ``|J|``: vector only, never the scalar magnitude.

**Reference**: commit ``8dbee35`` (BEM B/J), commit ``1aa66ee`` (FEM
B/J).

============================================================
## gmsh_viz — Visualization quality checklist (regression guard)
============================================================

**Symptom**: New GMSH outputs look noticeably worse than older
outputs from the same panel. Common regression patterns:

  - sphere appears black (volume mesh edges fill the screen)
  - field views are present but show flat color (no arrows)
  - high-order curving is lost; coil renders as faceted polygon
  - only B is exported; J is missing
  - the Open GMSH button enables but the .msh has no field data

**Root cause**: Each new export refactor drops one or two of the
quality knobs the previous version had. Without a checklist the
quality drifts down monotonically with every refactor.

**Rule — checklist for any new GMSH export from a panel**:

DO:

  - Write a **vector** view per physical quantity (B, J, A, ...).
    Use ``add_field("B", node_data, ncomp=3)``.
  - Per-vertex evaluate the CoefficientFunction directly:

    ```python
    arr = np.zeros((mesh.nv, 3))
    for vi, v in enumerate(mesh.vertices):
        val = cf(mesh(*v.point))
        arr[vi] = [float(getattr(x, "real", x)) for x in val[:3]]
    ```

    NodeData on the original vertices renders cleanly even on
    high-order curved meshes — no projection-induced artifacts.

  - Export BOTH the field-source view (J on the coil) AND the
    field-effect view (B everywhere). One view alone is half the
    story; the user always wants to correlate the two.

  - Restrict source-side views to the relevant material:
    ``add_field("J", node_J, ncomp=3, material="coil")`` so the
    arrows do not pollute the air domain.

  - Write a **companion ``.geo`` file** next to the .msh and have
    the panel's "Open GMSH" button open the .geo, NOT the .msh.
    The .geo applies the display options below automatically.

  - Apply these display options in the .geo:

    ```
    Mesh.NumSubEdges = 4;       // smooth high-order curving
    Mesh.Volumes = 0;           // hide air tetrahedra
    Mesh.VolumeFaces = 0;
    Mesh.VolumeEdges = 0;
    Mesh.SurfaceEdges = 0;      // hide internal triangle edges
    Mesh.SurfaceFaces = 0;
    View[0].VectorType = 4;     // 3D arrow style for B
    View[0].IntervalsType = 3;  // continuous color map
    View[1].VectorType = 4;     // 3D arrow style for J
    View[1].IntervalsType = 3;
    ```

  - Activate ``mesh.Curve(N)`` (N >= 2) on the source mesh BEFORE
    writing so the GMSH exporter can emit Tri6/Tri10 elements with
    correct mid-edge nodes. Linear-only export of a high-order .vol
    silently degrades the geometry visualization.

  - Always provide a **fallback** path that writes a mesh-only
    .msh when the field-write fails (see ``silent_except`` topic),
    so the Open GMSH button always enables.

DON'T:

  - **Do not** add a separate ``|B|`` / ``|J|`` scalar view. The
    vector view already encodes magnitude in its arrow color. The
    scalar duplicates information and clutters the View list.

  - **Do not** leave ``Mesh.Volumes = 1`` (the GMSH default). On a
    100k-element sphere mesh that draws as a black blob and hides
    every field arrow.

  - **Do not** write per-element ElementData when per-vertex
    NodeData is what GMSH wants for arrow rendering. Per-element
    arrows are constant within each cell and look chunky;
    per-vertex arrows interpolate smoothly.

  - **Do not** project a vector quantity to a scalar GridFunction
    just to use ``Set()``. ``H1(mesh, dim=3)`` is not a real API;
    ``VectorH1`` and ``H1(...)**3`` work but add nothing over
    direct CF evaluation, and ``Set(complex_cf)`` on a real space
    raises silently. Skip the projection entirely.

  - **Do not** ship just B without J (or vice versa). The user
    always wants both.

  - **Do not** open the raw .msh from the panel button. Open the
    companion .geo instead so the display options are applied.

**Reference**: commits ``d364796`` (FEM B/J vector + companion .geo
+ mesh-only fallback), ``1aa66ee`` (initial unified BEM/FEM viz),
``8dbee35`` (BEM dropping ``|B|``/``|J|`` scalars). Past quality
regressions traced from comparing ``calc_fem_kelvin.py`` Step 7
across these commits.

============================================================
## gmsh_arrow_size — Without ArrowSizeMin/Max the field is invisible
============================================================

**Symptom**: GMSH .geo opens, B view is selected, the user sees a
clean wireframe coil but **no arrows at all**. Or: arrows are
visible but tiny dots, indistinguishable from the background mesh.
The user reports "出力されるgmshの質が過去のものより劣化していた".

**Root cause**: GMSH's vector view defaults to ``ArrowSize* = 60``
**relative to the bounding box**, which for a 0.24 m air sphere
gives arrows that are physically a few mm long — invisible at the
default zoom. ``View.VectorType = 4`` alone is not enough; without
fixed pixel sizes the arrows do not render at a useful scale.

The v3.6.1 (commit ``5e0b69c``, 2026-03-23) version of
``calc_inductance.py post_process`` set:

```
View[0].ArrowSizeMin = 20;
View[0].ArrowSizeMax = 20;
View[1].ArrowSizeMin = 20;
View[1].ArrowSizeMax = 20;
```

That gave the consistently beautiful field arrows everyone praised.
A 2026-04-12 refactor of the same companion .geo writer dropped
these lines (left it at GMSH defaults), and the field views became
invisible without the user noticing — the comment in the new code
even said "the arrow size is left at GMSH's auto-scale" as if that
were a feature.

**Rule**: Every ``.geo`` companion file the panel writes MUST set:

```
View[N].VectorType = 4;          // 3D arrow style
View[N].IntervalsType = 3;       // continuous color map
View[N].ArrowSizeMin = 20;       // FIXED 20 px — do NOT auto-scale
View[N].ArrowSizeMax = 20;
```

for **every** vector view in the .msh. Without ArrowSizeMin/Max
the field is functionally invisible. There is no scenario where
GMSH's auto-scale produces a useful default for Radia output, so
the "leave it default" path is wrong by definition.

**Reference**: commit ``5e0b69c`` (original beautiful version),
commit ``d949d50`` (regression: dropped the arrow sizes), this
commit (restored).

============================================================
## subprocess_args — calc_*.py choices must match GUI combos
============================================================

**Symptom**: GUI combo has items ``[sibc, esim]``, but the calc
script's ``argparse`` declares ``choices=["esim", "dowell"]``. User
picks "sibc" in the GUI, the subprocess dies with
``argument --impedance-model: invalid choice: 'sibc'``.

**Root cause**: GUI labels and CLI labels drifted independently.
Often the GUI uses the more user-friendly term and the script keeps
the original internal name.

**Rule**: When you rename a combo item, EITHER:

  1. Add the new name as an explicit alias in the calc script's
     argparse ``choices`` and translate to the internal name in the
     ``__main__`` block (``imp == "sibc" -> imp = "dowell"``), OR
  2. Update both sides in the same commit and verify.

Either way, write a quick smoke test that constructs the command via
``_build_*_command`` and checks ``argparse.parse_args(cmd[1:])``
accepts every option.

**Reference**: commit ``8dbee35``, ``calc_inductance.py`` accepts
``sibc`` as alias for legacy ``dowell``.

============================================================
## subprocess_args — Hidden widgets feed argparse
============================================================

**Symptom**: GUI runs the FEM solver successfully but the result is
nonsense — sigma reports 2e6 even though user typed 5.8e7 in
the BEM tab.

**Root cause**: ``self.val("wp_sigma")`` returns the LAST text the
user typed in that widget, regardless of whether the widget is
currently visible. Mode switches don't reset values — they only
toggle visibility. So a value set under BEM persists across mode
switches and silently feeds the FEM call.

**Rule**: Either keep the widget VISIBLE in every mode that reads
it (preferred — see ``layout_unification``), or RESET its value to
the per-method default in ``_on_method_changed``. Never read a hidden
widget in ``_build_*_command``.

**Reference**: commit ``1aa66ee``, ``IHPanel._build_fem_command``
now reads only widgets that are visible in FEM mode.

============================================================
## cubit_jou — `subtract A from B keep_tool` id semantics
============================================================

**Symptom**: ``${kelvin_top = Id("volume")}`` after a subtract
returns the wrong volume id. Subsequent block / sideset commands
silently target the wrong volume, blocks come up empty, sidesets
fail with "no surface".

**Root cause**: Cubit 2025.3's ``subtract A from B keep_tool`` has
two behaviors depending on B's id:

  - **B is below the current max id**: B is MODIFIED in place. Its
    id stays the same. ``Id("volume")`` after the call returns the
    global max (some unrelated volume), not B.
  - **B IS the current max id**: B is REPLACED with a new id =
    max + 1. ``Id("volume")`` after the call returns the new id.

The two cases are silent — Cubit's "Updated volume(s): N" message
shows the modified id either way.

**Rule**: Capture volume ids RIGHT AFTER the operation that creates
them, and re-capture only when the operation creates a new id.
Single-tool subtract on a non-max body does NOT create a new id; do
not re-capture in that case. Multi-tool subtract on a max-id body
DOES create a new id; re-capture with ``${name = Id("volume")}``.

When in doubt, write a smoke ``test_*.jou`` that runs the script
through ``coreform_cubit -batch -nographics`` and verifies the block
membership.

**Reference**: commit ``ac64254``,
``validation_test/induction_heating/demoted_samples_legacy/ih_fem_kelvin_sample.jou``
header comment "subtract ... keep_tool id semantics".

============================================================
## cubit_jou — Sweep gap-face surface ids drift after imprint
============================================================

**Symptom**: ``sideset 1 add surface 1`` works in ``ih_bem_sample.jou``
(simple subtract sequence) but fails in
``validation_test/induction_heating/demoted_samples_legacy/ih_fem_kelvin_sample.jou``
(extra webcuts + subtracts) — surface 1 no longer exists after the
multiple imprint+merge cycles.

**Root cause**: Each ``imprint`` operation can renumber surface ids.
Hardcoded IDs are valid only for the EXACT operation order they were
captured under.

**Rule**: For sample .jou files that go through multiple imprint
cycles, identify the gap faces by **area + centroid** instead of by
hardcoded id:

```
group "coil_gaps" add surface in volume {coil_vid} with area < 0.0001
sideset 1 add surface in coil_gaps with y_coord > -0.001
sideset 1 name "source"
sideset 2 add surface in coil_gaps with y_coord < -0.001
sideset 2 name "sink"
```

Keep simpler scripts (single subtract, no webcut) on hardcoded ids
if they are verified to work — those are easier to read.

**Reference**: commit ``ac64254``,
``validation_test/induction_heating/demoted_samples_legacy/ih_fem_kelvin_sample.jou``
section 11 ``Sidesets``.

============================================================
## sample_jou — One sample per (panel, solver method) pair
============================================================

**Symptom**: A single ``ih_sample.jou`` was used for both BEM and FEM
modes of the IH panel, but BEM needs only a small air sphere
(surface integral equation, open BC) and FEM needs either a large
air sphere (Dirichlet truncation) or a Kelvin exterior domain (exact open BC).
Forcing one geometry to serve both made BEM slow and FEM inaccurate.

**Rule**: Ship one sample .jou per (panel, solver method) pair when
the methods need fundamentally different mesh topology. Naming
convention:

  ``{stem}_{method}_sample.jou``

Examples (current layout after the 2026-04-23 demotion):

  panels/samples/ih_bem_sample.jou        (BEM, shipped in wheel)
  panels/samples/ih_peec_bem_coarse.jou   (peec_bem canonical, golden)
  panels/samples/ih_fem_kelvin_skin_fine.jou
                                          (fem_coilmesh canonical, golden)
  panels/samples/ih_peec_inductance.jou   (peec_inductance canonical)

The no-Kelvin FEM baseline (``ih_fem_sample.jou``) and the small-mesh
FEM + Kelvin variant (``ih_fem_kelvin_sample.jou``) that used to live
next to these are now under
``validation_test/induction_heating/demoted_samples_legacy/`` — they were
not golden-tested and the Kelvin variant shipped with a
misleading "auto-Kelvin" comment.  Any NEW sample must stay aligned
with the (panel, solver method) = one-file rule and be golden-locked
before joining ``panels/samples/``.

**Reference**: commit ``ac64254`` (split), commit ``241090c``
(restored BEM-validated mesh density), 2026-04-23 demotion commit
(moved non-canonical samples out of shipped panel samples).

============================================================
## silent_action — Menu actions must produce visible feedback
============================================================

**Symptom**: User clicks "Verify Deploy" in the Solve menu. Nothing
happens — no dialog, no console output, no busy cursor. The user
reports "押しても何もおきない".

**Root cause**: The handler wrote its mtime/sha1 report to
``C:/radia_panel_log.txt`` via ``_panel_log()`` and **only** there.
The user has no reason to be tailing that file, so the action is
indistinguishable from a no-op.

**Rule**: Any user-triggered menu action MUST produce a visible
result. Acceptable channels:

  1. **QMessageBox** for ad-hoc reports / confirmations
  2. The panel **output text widget** (``self._output.appendPlainText``)
     for streaming subprocess output
  3. A new **window** for non-trivial reports

Writing to a debug log file is fine as a SECONDARY record but never
the only output. Test by clicking the button: if the screen does not
change, the handler is broken regardless of whether the file got
written.

**Reference**: commit ``d364796``, ``register_toolbar.py::
_verify_deploy``.

============================================================
## silent_except — Bare except hides field-export bugs
============================================================

**Symptom**: GMSH "Open" button stays disabled after a successful
FEM run. Result panel shows L, P, Q, time — everything looks normal.
Panel debug log says ``msh_file (value='') — Open GMSH disabled``.

**Root cause**: The GMSH export step inside ``calc_fem_kelvin.py``
was wrapped in ``try / except Exception as e: _log(f"GMSH_ERROR:{e}")``.
The actual exception was ``H1() got an unexpected keyword argument
'dim'`` (the ``H1(mesh, order=1, dim=3)`` API does not exist — vector
H1 is built via ``H1(...)**3`` or ``VectorH1``). The bare except
swallowed the AttributeError, ``_log`` printed only the short
``str(e)`` to stderr (which was already past the visible output
window), and ``gmsh_file`` stayed at its initial empty string.

The result dict shipped ``msh_file=""``, the panel saw an empty
path, and disabled the Open GMSH button — looking exactly like a
"no field exported, success ignored" silent failure.

**Rules**:

  1. Never log just ``str(e)`` from a bare except. Log
     ``type(e).__name__`` AND the **last few lines of the
     traceback** so the actual offending API call is visible:

     ```python
     except Exception as e:
         import traceback
         tb = traceback.format_exc()
         _log(f"GMSH_ERROR:{type(e).__name__}: {e}")
         for line in tb.splitlines()[-4:]:
             _log(f"GMSH_ERROR:  {line}")
     ```

  2. **Always provide a fallback** that still produces SOMETHING
     openable. Even if the field-write fails, write a mesh-only
     ``.msh`` so the Open GMSH button enables and the user gets
     the geometry:

     ```python
     except Exception:
         post_min = GmshPostExport(mesh, boundary=False)
         post_min.write_mesh(msh_output)
         gmsh_file = msh_output
     ```

  3. Don't use ``H1(mesh, dim=N)`` — that signature does not exist
     in NGSolve. For per-vertex vector evaluation in a GMSH export,
     skip the GridFunction projection entirely and evaluate the
     CoefficientFunction directly at each ``mesh.vertices[i].point``.
     The (n_v, 3) numpy array goes straight into the GMSH NodeData
     section.

**Reference**: commit ``d364796``, ``calc_fem_kelvin.py`` Step 7.

============================================================
## result_keys — Subprocess result dict keys are an API contract
============================================================

**Symptom**: One panel script returns ``"msh_file"`` in its result
dict, another returns ``"msh_output"``, a third returns
``"gmsh_file"``. The GUI's ``_on_finished`` has to test for all of
them, and any new script that picks a different name silently
breaks the Open GMSH button.

**Root cause**: Each calc_*.py was developed independently and
each author picked a slightly different key name for the same
concept ("path of the GMSH file the user should open").

**Rule**: Subprocess result dicts have a **stable, documented set
of keys** that the GUI relies on. The current contract for the IH
panel:

  ``msh_output`` or ``msh_file``  -- raw .msh path the user can open
                                     (Open GMSH button uses either)
  ``field_gmsh_file``              -- preferred .geo path that Merges
                                     mesh + field views with display
                                     options applied (overrides
                                     msh_output when present)
  ``inductance_H``                 -- headline result for L
  ``coupled_*``                    -- coupled BEM headline results
  ``error``                        -- single string, presence
                                     suppresses success display

When you add a new key, update the gui_base.py ``_on_finished``
handler in the same commit. When you rename a key, grep for it
across ``radia_gui_base.py`` AND every ``radia_*.py`` panel before
shipping.

**Reference**: ``radia_gui_base.py::_on_finished`` lines 425-470,
which has the canonical "what the GUI looks for" logic.

============================================================
## regression_blast_radius — One edit breaks an unrelated place
============================================================

**Symptom**: Fixing the FEM panel's Open GMSH button suddenly
breaks the BEM panel. Adding a return_vertex_map kwarg to
``_extract_surface_mesh_filtered`` makes BEM die with
``TypeError: int() argument must be a string ... not PointId``
on a code path the FEM fix did not touch. The user reports
"どっかをいじるとどっか別のところが壊れる".

**Root cause**: The Radia panels share a lot of utility code via
``calc_common.py`` / ``calc_heating_bem.py`` / the BEM solver. A
"local" change to one panel that touches a shared helper has a
**blast radius** that the immediate test doesn't catch:

  - Adding a kwarg with a new return type breaks every existing
    caller that did not opt in.
  - Renaming a result-dict key breaks the GUI's ``_on_finished``
    consumer regardless of which panel script wrote it.
  - Changing a calc_*.py argparse default breaks any sample run
    that relied on the old default being passed implicitly.
  - Tightening a CoefficientFunction signature breaks every
    GridFunction.Set() call in the same module.

The PointId crash above came from a one-line construction:

```python
new_to_old = {int(new_id) - 1: int(old_nr)
               for old_nr, new_id in old_to_new.items()}
```

``int(PointId)`` raises TypeError on netgen's PointId class.
Discovered only when the user ran BEM (not FEM, the panel I was
editing).

**Rules**:

  1. **Run BOTH panels** after touching any shared helper, even
     when "the change is only for FEM". If the helper lives in
     ``calc_heating_bem.py`` or ``calc_common.py``, the smoke
     test must cover at least one BEM run AND one FEM run.

  2. **Optional kwargs** that change return type are dangerous —
     prefer adding a separate function or making the new return
     a tuple ``(old_result, extra)`` with a clear opt-in flag,
     and write a regression test that exercises the OLD signature.

  3. **Cast carefully** — ``int(x)`` works for built-in numerics
     but not for opaque library wrappers like netgen's PointId.
     When in doubt use the explicit attribute (``x.nr`` for
     PointId) and document why in a comment.

  4. **Result dict keys** are an API contract (see ``result_keys``
     topic). Renaming requires grep across every consumer.

  5. After any non-trivial edit, **list every file that imports
     the changed function** and run the smoke test for each
     consumer panel. ``grep -rln "_extract_surface_mesh_filtered"
     src/radia/`` is faster than discovering the regression in
     production.

**Reference**: commits ``1aa66ee`` (introduced
``return_vertex_map``), ``d364796`` (broke BEM via PointId int
cast), this commit (restored ``.nr`` access + added rule).

============================================================
## panel_qt_testing — retired; notebook tests are the active gate
============================================================

This topic is historical.  The PySide6 desktop panel route has been retired in
favor of Jupyter notebook workbenches.  Current panel behavior is checked by:

```
python -m pytest validation_test/panels/test_notebook_workbench.py -q
```

The active invariants are:

- notebooks do not import PySide6 / PyQt
- `DesignSpec(...)` cells are the canonical initial-value store
- JSON files are run artifacts, not presets
- `CommandWorkbench.run_local()` writes `radia_result.v2`
- `Workbench.build_command()` stays aligned with the target `calc_*.py`

Normal Radia Python should not install PySide6.  Coreform Cubit's embedded
PySide6 is protected and must not be removed.

============================================================
## learn_edition_cap — Cubit Learn Edition 50k limit is harmless
============================================================

**Symptom**: Running a sample .jou through Cubit Learn Edition prints
``ERROR: Coreform Cubit - Learn Edition restricts export to models
with less than 50k elements.`` for a 147k-element model. New
contributors think the export failed and start coarsening the mesh.

**Root cause**: The cap applies to Cubit's own ``export gmsh`` /
``export vtk`` / ``export exo`` exporters. The Radia in-tree
``export netgen`` plugin BYPASSES the cap and writes the .vol
regardless of the warning. Both the ERROR line and the
``Exported Netgen Vol (order N): ...`` line appear in the same Cubit
run, so the export succeeded.

**Rule**: Tune sample .jou meshes to whatever density the physics
needs (typically the BEM-validated reference density). Never coarsen
to "fit under 50k". The deployed LAB / 100号機 / mdx machines all
run Cubit Pro and never see the warning.

**Reference**: ``CLAUDE.md`` § "Cubit Learn Edition 50k Element Cap",
commit ``334d84f``.

============================================================
## wheel_guard — Mouse wheel silently changes combo/spin values
============================================================

**Symptom**: The user scrolls the panel with the mouse wheel to reach
a field further down. As the cursor passes OVER a combo box (e.g. the
one just below "Workpiece impedance model") or a spin box, that widget
**eats the wheel event and changes its own selection** instead of
scrolling the page. The user does not notice — they ran the solve
with ``esim`` selected when they had carefully picked ``sibc`` a
moment earlier. The wrong value is invisible because the combo looks
identical; only the result comes out wrong.

**Root cause**: Qt's default ``QComboBox`` / ``QAbstractSpinBox``
``wheelEvent`` steps the value whenever the pointer is over the widget
(no focus required). In a tall scrolling form this is almost always
the WRONG behaviour — the user means to scroll the form, not to nudge
a parameter. Because the change is a single step (one combo item, one
spin increment) it produces a plausible-looking value, so it is the
worst kind of silent UI bug: the same class as a silent fallback (a
plausible wrong number the user cannot audit — see "No Fallbacks").

**Rule**: Every value-bearing widget in a scrollable panel MUST
swallow wheel events. Use the ``_NoWheel*`` subclasses in
``radia_gui_base.py``:

```python
class _NoWheelComboBox(QComboBox):
    def wheelEvent(self, event):
        event.ignore()          # let the QScrollArea scroll instead

class _NoWheelSpinBox(QSpinBox): ...
class _NoWheelDoubleSpinBox(QDoubleSpinBox): ...
```

``ModePanel.add_combo`` builds ``_NoWheelComboBox`` and ``add_spin``
builds ``_NoWheelSpinBox``, so every panel field declared through the
base class gets the guard for free. Any panel that constructs a raw
``QComboBox`` / ``QSpinBox`` / ``QDoubleSpinBox`` directly in a temporary
legacy adapter MUST import the ``_NoWheel*`` alias instead of the bare Qt
class:

```python
from radia_gui_base import (_NoWheelComboBox as QComboBox,
                            _NoWheelSpinBox as QSpinBox,
                            _NoWheelDoubleSpinBox as QDoubleSpinBox)
```

**Test**: ``tests/panels/test_combo_wheel_guard.py`` synthesises a
``QWheelEvent`` over every combo / spin found via ``findChildren`` in
each panel window and asserts the value did NOT change. This is a
FUNCTIONAL test (it dispatches the real event), not a string grep —
a panel that forgets the guard fails it.

**Reference**: ``radia_gui_base.py`` ``_NoWheel*`` classes +
``add_combo`` / ``add_spin``; ``CLAUDE.md`` "Panel Layout Policy"
(2026-05-29); the ``panel-wheel-guard`` skill.

============================================================
## panel_layout — 10pt font + vertical scrollbar, never compress rows
============================================================

**Symptom**: On a short window (laptop screen, or the user dragged the
panel small) the parameter rows get squeezed vertically until the
field text is **clipped** — the entered numbers are half-cut and
unconfirmable. The user reports "行が細くなりすぎてて、文字潰れ"
(rows too thin, text crushed).

**Root cause**: A bare ``QFormLayout`` / ``QVBoxLayout`` packed
straight into the window will compress its rows below their natural
``sizeHint`` height when the window is shorter than the content. Qt
shrinks every row proportionally rather than introducing a scrollbar,
so a 20-row form in a 400px window clips every row.

**Rule (Panel Layout Policy, 2026-05-29)**:

  1. Base panel font is **10pt** — ``PANEL_BASE_FONT_POINT_SIZE = 10``
     in ``radia_gui_base.py``, applied to the panel widget. Do NOT
     shrink below 10pt to "make it fit".
  2. The parameter form is wrapped in a **``QScrollArea`` with
     ``setWidgetResizable(True)``**. When the window is too short the
     area shows a **vertical scrollbar** and every row keeps its
     natural (uncompressed) ``sizeHint`` height. NEVER let rows
     compress.
  3. The horizontal scrollbar is disabled
     (``setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)``) and the
     frame is ``QFrame.NoFrame`` so the scroll area is visually
     transparent.

```python
scroll = QScrollArea()
scroll.setWidgetResizable(True)              # rows keep sizeHint height
scroll.setFrameShape(QFrame.NoFrame)
scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
scroll.setWidget(form_container)
```

This pairs with ``wheel_guard``: once the form scrolls, the
``_NoWheel*`` widgets ensure the wheel scrolls the area instead of
nudging whatever combo the cursor happens to be over.

**Test**: ``tests/panels/test_panel_output_health.py``
``test_panel_base_font_is_10pt`` (asserts the constant) and
``test_form_is_in_resizable_scrollarea`` (asserts a
``widgetResizable`` ``QScrollArea`` exists in the window).

**Reference**: ``radia_gui_base.py`` ``PANEL_BASE_FONT_POINT_SIZE`` +
``_build_ui`` QScrollArea wrap; ``CLAUDE.md`` "Panel Layout Policy:
10pt + Vertical Scrollbar, Never Compress" (2026-05-29).

============================================================
## result_output — OUTPUT/JSON must surface ne / DoF / time + integral quantities
============================================================

**Symptom**: A solve finishes, the result panel shows the headline
number (inductance, say) but **no element count, no DoF, no
compute-time breakdown, no heat, no temperature**. The user cannot
tell how big the model was or how long each phase took. For an IH run
the heat (P) and the workpiece temperature — the whole point of the
analysis — are missing from both the OUTPUT text and the saved JSON.

**Root cause (two distinct bugs)**:

  1. ``_on_finished``'s summary block checked a FIXED cascade of key
     names (``n_dofs`` / ``t_solve`` / ``P_total_W``) while the
     calc_*.py scripts actually emit ``wp_ndof`` / ``t_bem_solve_s`` /
     ``P_wp_W``. The names never matched, so the summary rendered
     nothing — the data was IN the result dict, just never displayed.
  2. The thermal script reported a single peak temperature, not the
     mean / max / min triple the user needs to judge a heat-treat
     result.

**Rule (Result Output Policy, 2026-05-29)**:

  1. The OUTPUT text AND the ``--output`` JSON MUST record, for every
     analysis: **element count, DoF, and a per-phase compute-time
     breakdown**. ``calc_common.calc_main`` already writes the full
     result dict to the JSON; the GUI must surface the same keys.
  2. ``_append_standard_summary(result)`` (in ``radia_gui_base.py``)
     renders these **generically from whatever keys are present** — it
     does NOT hard-code a key cascade. It surfaces element / triangle
     counts (``ne`` / ``wp_mesh_n_tris``), DoF (``ndof`` / ``wp_ndof``),
     every ``t_*_s`` timing, heat ``P_*_W``, and temperature.
  3. For IH specifically, **heat and temperature are mandatory**.
     Temperature is reported as **mean (volume-averaged) / max / min**,
     not a single peak — the calc emits ``T_mean_C`` (=
     ``Integrate(gfT, wp_mesh) / volume``), ``T_max_C``, ``T_min_C``.
  4. When you add a new integral quantity to a calc script it appears
     in the OUTPUT automatically IF you name it with a recognised
     convention (``t_*_s`` for time, ``P_*_W`` for power, ``*_C`` for
     temperature, ``ndof`` / ``ne`` for sizes). Extend
     ``_append_standard_summary`` rather than re-implementing the
     summary per panel.

**Why generic, not a key cascade**: a fixed list of key names is
exactly the bug above — every new calc script that picks a slightly
different name silently falls out of the summary (same failure mode
as ``result_keys`` for the Open GMSH button). Rendering from the keys
that are present means a new ``t_assembly_s`` or ``P_coil_W`` shows up
with zero GUI changes.

**Test**: ``tests/panels/test_panel_output_health.py``
``test_output_summary_surfaces_ne_dof_time_heat`` (BEM-shaped mock:
asserts Elements / DoF / Heat (P) / Compute time + the actual key name
``t_bem_solve_s`` all appear) and
``test_output_summary_temperature_mean_max_min`` (thermal mock:
asserts "Temperature [C]" with mean / max / min values).

**Reference**: ``radia_gui_base.py::_append_standard_summary`` +
``_on_finished``; ``calc_heat.py`` ``T_mean_C`` / ``T_max_C`` /
``T_min_C``; ``CLAUDE.md`` "Result Output Policy: ne / DoF / time +
analysis integral quantities" (2026-05-29).
"""


_TOPICS = (
    "combo_state",
    "mode_switch",
    "layout_unification",
    "gmsh_viz",
    "gmsh_arrow_size",
    "subprocess_args",
    "cubit_jou",
    "sample_jou",
    "silent_action",
    "silent_except",
    "result_keys",
    "regression_blast_radius",
    "panel_qt_testing",
    "learn_edition_cap",
    "wheel_guard",
    "panel_layout",
    "result_output",
)


def get_panel_gui_pitfalls(topic: str = "") -> str:
    """Return the panel GUI pitfalls knowledge.

    Args:
        topic: Empty for the full document, or one of the entries in
               ``_TOPICS`` above.

    Multiple sections share the same ``gmsh_viz`` keyword (vector-only
    + companion .geo); requesting that topic returns BOTH sections so
    the reader gets the complete pitfall together.
    """
    if not topic:
        return PANEL_GUI_PITFALLS

    if topic not in _TOPICS:
        return (f"Unknown topic: {topic!r}. Available topics:\n"
                f"  {', '.join(_TOPICS)}\n\n"
                f"Pass empty string for the full document.")

    # Build a list of (topic_keyword, char_offset) for every section
    # header in the document. The header line format is exactly
    # ``## <topic> — <title>`` (one space, em-dash). Then slice from
    # the requested topic's previous "## ===" delimiter up to the
    # NEXT topic's "## ===" delimiter (or end of document).
    headers = []
    pos = 0
    while True:
        next_pos = PANEL_GUI_PITFALLS.find("\n## ", pos)
        if next_pos < 0:
            break
        line_end = PANEL_GUI_PITFALLS.find("\n", next_pos + 1)
        line = PANEL_GUI_PITFALLS[next_pos + 1:line_end]
        for t in _TOPICS:
            if line.startswith(f"## {t} "):
                headers.append((t, next_pos + 1))
                break
        pos = next_pos + 1

    # Find requested + next-distinct-topic positions.
    req_starts = [off for kw, off in headers if kw == topic]
    if not req_starts:
        return f"Topic {topic!r} declared but not found in document."
    section_start = PANEL_GUI_PITFALLS.rfind(
        "## ===", 0, req_starts[0])
    if section_start < 0:
        section_start = req_starts[0]

    # End = start of the next section whose topic differs from the
    # requested one. This keeps both ``gmsh_viz`` sections together.
    section_end = len(PANEL_GUI_PITFALLS)
    last_req = req_starts[-1]
    for kw, off in headers:
        if kw != topic and off > last_req:
            # Slice ends at the delimiter line just above this header.
            delim = PANEL_GUI_PITFALLS.rfind("## ===", 0, off)
            section_end = delim if delim > 0 else off
            break

    return PANEL_GUI_PITFALLS[section_start:section_end].rstrip() + "\n"
