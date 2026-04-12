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
cubit_jou) so they can be referenced from prompts.
"""

PANEL_GUI_PITFALLS = """\
# Radia GUI / Cubit Panel Pitfalls

This document lists pitfalls we have hit (and fixed) while developing
the Radia panels and the calc_*.py subprocess scripts. Each section
ends with a **rule** to prevent the same bug from coming back.

Topics: combo_state, mode_switch, gmsh_viz, subprocess_args, cubit_jou,
        sample_jou, layout_unification, learn_edition_cap

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

**Reference**: commit ``1aa66ee``, ``radia_ih.py::IHPanel._build_ui``
("Layout (per user request 2026-04-12)" comment block).

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

**Reference**: commit ``ac64254``, ``ih_fem_kelvin_sample.jou``
header comment "subtract ... keep_tool id semantics".

============================================================
## cubit_jou — Sweep gap-face surface ids drift after imprint
============================================================

**Symptom**: ``sideset 1 add surface 1`` works in ``ih_bem_sample.jou``
(simple subtract sequence) but fails in ``ih_fem_kelvin_sample.jou``
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

**Reference**: commit ``ac64254``, ``ih_fem_kelvin_sample.jou``
section 11 ``Sidesets``.

============================================================
## sample_jou — One sample per (panel, solver method) pair
============================================================

**Symptom**: A single ``ih_sample.jou`` was used for both BEM and FEM
modes of the IH panel, but BEM needs only a small air sphere
(surface integral equation, open BC) and FEM needs either a large
air sphere (Dirichlet truncation) or a Kelvin shell (exact open BC).
Forcing one geometry to serve both made BEM slow and FEM inaccurate.

**Rule**: Ship one sample .jou per (panel, solver method) pair when
the methods need fundamentally different mesh topology. Naming
convention:

  ``{stem}_{method}_sample.jou``

Examples:

  ``ih_bem_sample.jou``         BEM, small air, surface mesh
  ``ih_fem_sample.jou``         FEM, large air sphere, Dirichlet
  ``ih_fem_kelvin_sample.jou``  FEM + Kelvin shell, exact open BC

Update ``CONVENTIONS.md`` and ``panel_conventions_knowledge.py`` so
the registry knows about all three.

**Reference**: commit ``ac64254`` (split), commit ``241090c``
(restored BEM-validated mesh density).

============================================================
## learn_edition_cap — Cubit Learn Edition 50k limit is harmless
============================================================

**Symptom**: Running a sample .jou through Cubit Learn Edition prints
``ERROR: Coreform Cubit - Learn Edition restricts export to models
with less than 50k elements.`` for a 147k-element model. New
contributors think the export failed and start coarsening the mesh.

**Root cause**: The cap applies to Cubit's own ``export gmsh`` /
``export vtk`` / ``export exo`` exporters. The Radia in-tree
``radia_export netgen`` plugin BYPASSES the cap and writes the .vol
regardless of the warning. Both the ERROR line and the
``Exported Netgen Vol (order N): ...`` line appear in the same Cubit
run, so the export succeeded.

**Rule**: Tune sample .jou meshes to whatever density the physics
needs (typically the BEM-validated reference density). Never coarsen
to "fit under 50k". The deployed LAB / 100号機 / mdx machines all
run Cubit Pro and never see the warning.

**Reference**: ``CLAUDE.md`` § "Cubit Learn Edition 50k Element Cap",
commit ``334d84f``.
"""


_TOPICS = (
    "combo_state",
    "mode_switch",
    "layout_unification",
    "gmsh_viz",
    "subprocess_args",
    "cubit_jou",
    "sample_jou",
    "learn_edition_cap",
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
