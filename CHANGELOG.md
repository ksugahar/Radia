# Changelog

All notable changes to the `radia` package.  Format: each release lists
**what shipped** + **why** in compact form.  Packaged wheels on PyPI.

## 4.95.0 — Multipole-moment MMM matrix-free and HACApK acceleration

Released 2026-06-26.

- **Multipole-moment MMM method 1**: replaced the dense BiCGSTAB linear
  step with matrix-free moment matvecs and element-wise block Jacobi for
  pure hexahedral and mixed hex/wedge/pyramid 5/6-DOF surface-charge
  systems.  The path now uses the same moment blocks on demand instead of
  building a full dense matrix.
- **Multipole-moment MMM method 2**: the nonlinear HACApK path now reuses a
  chi-free geometry H-matrix across Picard iterations and applies the
  current per-element susceptibility as `Lx + diag(chi) Kx`, rebuilding
  only the local block-Jacobi preconditioner and RHS each outer iteration.
  The default Krylov solver remains BiCGSTAB; restarted GMRES is exposed as
  an explicit comparison path.
- **Fail-loud cleanup**: removed the inexact BiCGSTAB and two-level coarse
  preconditioner experiments from the runtime API.  Passing their old
  `SolverConfig` keys now raises instead of silently selecting a different
  path.  The failed branches are recorded in
  `docs/multipole_moment_mmm/MEMORY.md`.

## 4.89.1 — radia-ih hardening + stream_function Tikhonov/Pareto

Released 2026-06-02.

- **radia-ih**: removed the orphan `calc_heating.py` (legacy 2D-axisym
  FEM-ESIM that auto-generated its own mesh, violating the .vol-only
  Layer-4 contract; unreachable from the panel).  Verified
  `calc_fem_coilmesh --impedance-model esim` (nonlinear Karl) runs +
  converges and dropped its stale "WIP, raises" labels (calc help +
  panel combo).  Added Cubit-free axisym + esim regression goldens.
- **stream_function**: `RegularizedTSVD` folds Tikhonov (alpha)
  regularization; new Pareto-front demos (geometry / NSGA-II / cylinder)
  for the (homogeneity, peak-J) trade-off.
- Docs / radia-mcp knowledge reframed from the removed auto-mesh
  "2-Stage" path to the current EM -> q_surf -> Thermal flow.

## 4.89.0 — Loop-free solution by default (Helmholtz-Hodge loop removal)

Released 2026-06-02.

The MSC/MMM solvers now return a **loop-free physical magnetization** by
default.  The discrete operator `N` has an exact null space of circulating
surface charges ("loops" = the cycle space of the element-adjacency graph);
these are non-physical (they produce no field) but at high permeability they
dominate the raw solved `sigma` (up to ~99% of its norm).  After the solve
converges, a Helmholtz-Hodge projection removes the loop component:
`c` solves `(L^T L) c = L^T sigma` (CG; `L` = topological cycle basis), then
`sigma -= L c`.

- **Default ON** (`rad.SetLoopProjection(True)` is the default; pass `False`
  to keep the raw loop-containing `sigma`).
- **Field-transparent.** `N L = 0`, so `N*sigma` -- the field -- is unchanged
  (`rad.Fld` differs by `~1e-15`).  Verified: 102/102 core field + golden tests
  unchanged with the default on.  Only the (non-physical) circulating part of
  the magnetization distribution changes.
- **Cheap.** `L^T L` is the geometric, sparse, `mu_r`-independent, well-
  conditioned (`cond ~ 1`) loop Gram matrix; the CG converges in a handful of
  iters with no measurable slowdown.
- **All solver paths, consistent.** Method 0 (LU), 1 (dense BiCGSTAB), and 2
  (HACApK) all return the same loop-free `sigma` (the non-HACApK paths use a
  standalone pure-sparse cycle projection -- no H-matrix needed).
- **Nonlinear-safe.** `L = ker(N)` is `chi`-independent; the projection is
  applied once after the nonlinear iteration converges, so `chi(H)` is driven
  by the true loop-included `sigma` and only the final answer is loop-free.
  Auto-skipped when the loop-star / loop-deflated gauge is active.

Diagnostics: `rad.GetLoopProjStats()`.  Test: `tests/test_loop_projection.py`.
Docs: `docs/solver/MSC_NULLSPACE_DEFLATION.md` Section 8.5.  Also lands (opt-in,
off by default) the loop-deflated block-Jacobi gauge and the A_SS H-ILU
preconditioner as experimental reference paths.

## 4.88.0 — Cubit plugin command verb `radia_export` -> `export`

Released 2026-06-01.

The Cubit mesh-export commands now extend Cubit's native `export` verb
(`export netgen / gmsh / vtk / femeem / meg`, plus `export jmag_nastran`
for the JMAG-targeted Nastran BDF) instead of the old `radia_export`
verb.  The radia-side panels (`radia_export_menu`) emit the new commands;
the plugin binary lives solely in the cubit-mesh-export package (Tier-2),
so the radia wheel no longer bundles the `.ccm` / `.pyd`.

**Breaking**: existing `.jou` scripts calling `radia_export ...` must be
updated to `export ...` (the old verb is removed -> Cubit reports
`Unrecognized Keyword: 'radia_export'`).  See cubit-mesh-export 0.11.0
for the full command map.

## 4.87.0 — feat(ih): thermal targets the workpiece solid only (+ 2-file qsurf input)

Released 2026-06-01.

### feat(ih): thermal analysis targets the WORKPIECE SOLID only

`calc_heat.py` and `calc_heat_axisym.py` now REQUIRE a workpiece-only
volume mesh for the radia-ih Thermal step and reject anything else with
a clear, actionable error:

- **multi-material** (`>1` volume region, e.g. a coil+workpiece EM mesh)
  → "thermal analysis targets the WORKPIECE ONLY … has N volume
  materials […]"
- **surface-only** (`ne == 0`, a SIBC / PEEC surface mesh) → "needs a
  VOLUME mesh of the workpiece solid"
- **wrong dimensionality** (3D calc fed a 2D mesh, or vice versa) →
  points at the correct script

**Why**: the thermal mesh is, by design, a SEPARATE workpiece solid (the
EM mesh carries the workpiece as a SIBC-faced hole, WP-HOLE policy).
The solve previously ran on the WHOLE input mesh with no domain
restriction, so a coil+workpiece mesh silently diffused heat through the
coil as if it were steel (keiko 2026-05-31).  This is strict-in-what-we-
accept per the No-Fallback policy.  The dead `--material-label` flag
(declared but never used to restrict the domain) is removed.

The HeatPanel also warns on Browse (wp_vol tooltip) when the selected
.vol carries more than one volume material, so the mismatch is visible
before a Run round-trip.

Tests: `tests/panels/test_heat_workpiece_only.py` (multi-material
rejected, single-material accepted).  Also fixes a `from __future__`
import-ordering SyntaxError in `tests/panels/fixtures/generate_heat_cylinder.py`.

### feat(ih): thermal spatial source = qsurf .sol + workpiece .vol (em .vol auto-fills)

The HeatPanel spatial-q_surf input is now effectively two files — the
surface heat-density `qsurf .sol` and the workpiece volume `.vol`.
Picking the `.sol` auto-fills the EM `.vol` field from its companion
`<stem>_fem.vol` (the pair `calc_fem_kelvin.py` writes), shown in the
field and overridable.  A user-supplied EM `.vol` is never clobbered; a
`.sol` with no derivable companion clears a stale auto value so the user
browses the correct mesh.  The calc script still receives `--em-vol`
explicitly — the panel auto-fill is a visible UX convenience, not a
silent fallback (calc-side No-Fallback intact).

## 4.86.0 — feat: loop-star MSC gauge KEEPS the loops (field-exact, tree-cotree)

Released 2026-05-31.

`rad.SetLoopStarGauge(True)` (the HACApK/MSC tree-cotree loop-star solver)
now **keeps** the loop (null-space) content of the magnetic-surface-charge
solution instead of dropping it, so the external field matches the direct
LU / plain-BiCGSTAB / FEM solution to ~1e-10 at every `mu_r`.

**Why**: the old behaviour solved only the reduced star system
(`sigma = S y_S`) and discarded the loop part.  On the C-type electromagnet
the topological loops are not exactly externally field-silent *and* the
sparse star/loop bases are not exactly mutually orthogonal, so the external
field came out ~0.5% off (`dBz/Bz` = 5.4e-3 at `mu_r`=2) — a real method
error, not noise.

**Fix** (`RadHACApKMSCManager::SolveLoopStar`): after the reduced star solve,
recover the full solution by a few **block Gauss-Seidel** sweeps — (i) star
correction via the K-dense reduced solve, (ii) loop correction
`A_LL y_L = L^T r` with `A_LL = L^T diag(inv_chi) L` (CG) — each sweep
targeting the true residual `b - A sigma`.  Converges to the direct solution
regardless of star/loop basis orthogonality; the loops are kept.

**Verified field-exact** (`Bz` vs plain BiCGSTAB): C-type 6^3 `dBz/Bz`
4.7e-11…7.3e-10 across `mu_r` 2…1000 (was 5.4e-3…1.3e-5); cube 4^3 ~1e-11
(GS converges in 1 sweep); antisym IMA `mu_r`=1e5 `'+x'` 1.9e-4 / `'+x-z'`
7.4e-5 (both PASS, the high-`mu_r` shielding floor).

**Scope**: the star correction uses the K-dense `A_SS` LU (caps ~15^3 /
8 GB); field-exact keep-loops is validated in that regime.  Linear
(uniform-`chi`) target.  New diagnostic `rad.GetKeepLoopStats()`.  Larger-`N`
needs a scalable `A_SS` preconditioner (ILU/H-LU), the next increment.

## 4.85.3 — fix: linear-SIBC Z_s double-counted mu_r (extra sqrt(mu_r))

Released 2026-05-31.

`calc_inductance.py`'s Linear-SIBC workpiece branch multiplied
`(1+j)*rho/delta_wp` by an extra `math.sqrt(mu_r)`.  But `delta_wp`
(`EMMaterial.skin_depth`) already includes mu_r, so the surface impedance
came out proportional to mu_r instead of the standard Leontovich
sqrt(mu_r) -- i.e. sqrt(mu_r)x too large for magnetic (mu_r != 1)
workpieces.  Found via a radia-ih <-> COMSOL cross-validation (COMSOL
material sigma=5e5, mu_r=1000 gave radia Z_s=1.312 vs the correct 0.04149,
a factor sqrt(1000)~31.6; P_wp ~1/8).

The fix drops the extra `* math.sqrt(args.mu_r)`; the result now matches
`analytical_formulas.planar_surface_impedance` exactly and COMSOL.  Scope:
only mu_r != 1 (magnetic) workpieces in Linear-SIBC mode on the BEM-A /
PEEC weak-coupling inductance path.  Cu/Al (mu_r=1) were unaffected
(sqrt(1)=1); the ESIM and FEM paths use a different formula and were
unaffected.

Also in this release (panel maintenance, 2026-05-31): radia_em QSettings
restore-order fix (the launcher's .vol now wins over a stale saved
wp_vol), radia_motor Lamination "global" mode exposes --em-table, the
panel-cli-diff checker is generator-aware, and stale radia_accel /
radia_heat references were cleaned (those standalone panels were folded
into radia_em / radia_ih).

## 4.78.1 — vol/sol viewer: register via winreg (quote-mangling fix)

Released 2026-05-25.

Fixes `radia-vol-viewer --register` silently failing to install the
ftype handler on Windows.  The previous implementation passed an
`ftype` argument with embedded quotes through `subprocess.run(["cmd",
"/c", "ftype", 'Radia.VolViewer="C:\\...\\app.exe" "%1"'])`; the argv
→ Windows-command-line round-trip mangled the quotes so cmd.exe
received `ftype Radia.VolViewer=\C:\...\app.exe\ \%1\` and discarded
the registration.  Symptom: `--register` prints "Done", `.vol=Radia
.VolViewer` is set, but `cmd /c "ftype Radia.VolViewer"` returns
"file type not found" so double-click does nothing.

The fix replaces both `register_associations()` and
`unregister_associations()` with direct `winreg` writes to
`HKLM\SOFTWARE\Classes` (same scope as the original `assoc`/`ftype`,
no cmd.exe involvement, no quoting ambiguity).  `winreg` is stdlib;
no new dependency.

Discovered 2026-05-24 during the Stage 2 deploy to 100号機 + mdx,
worked around at the time with a manual `reg add` per-machine.

## 4.76.1 — .sol viewer: colour-mapped field (scalfunction/vecfunction)

Released 2026-05-24.

Fixes the v4.76.0 .sol display path which loaded the GridFunction and
called `Draw(gfu, mesh, name)` but left the togl viewer in
`visoptions.scalfunction = 'none'` -- the mesh appeared but with no
field colour overlay (the surfaces stayed uniform-green).  The fix
adds the missing Tcl bindings inside `view_sol._show()`:

```python
if is_vector:
    g.win.tk.eval(f'set visoptions.scalfunction "{name}:0"')
    g.win.tk.eval(f'set visoptions.vecfunction "{name}"')
    g.win.tk.eval('set visoptions.evaluate "abs"')
else:
    g.win.tk.eval(f'set visoptions.scalfunction "{name}"')
    g.win.tk.eval('set visoptions.vecfunction "none"')
g.win.tk.eval('set visoptions.showsurfacesolution 1')
g.win.tk.eval('Ng_Vis_Set parameters')
g.win.tk.eval('redraw')
```

The `_infer_fes()` helper now returns `(fes, label, is_vector)` so
the show callback can pick the right Tcl bindings.  HDiv / HCurl /
H1 dim=3 are flagged vector (colour = |F|, with arrow overlay), H1
scalar is flagged scalar (colour = value directly).

Probed via inspection of `visoptions.*` Tcl globals -- both
`scalfunction` and `vecfunction` start at `'none'` until set
explicitly, even after `Draw()` registers the scene.

## 4.76.0 — .vol / .sol double-click viewer (`radia-vol-viewer`)

Released 2026-05-24.

User-side replacement for the upstream NGSolve PR
[#242](https://github.com/NGSolve/netgen/pull/242) that adds
auto-load-on-launch for `.vol` files.  PR #242 was closed by the
author (head repo deleted), NOT rejected by maintainers; the upstream
behavior remains "open Netgen GUI with a blank viewport when invoked
on a `.vol`".  This release ships the workaround inside the wheel:

* New module `radia.tools.vol_sol_viewer` and two console scripts:
  - `radia-vol-viewer` — CLI (use for `--register` / `--unregister`
    and standalone invocations).
  - `radia-vol-viewer-gui` — `gui_scripts` variant launched via
    pythonw.exe (no console flash on double-click).  This is what
    `--register` writes into the Windows `.vol` / `.sol` registry
    associations.

* The viewer replicates upstream PR #242's `_vol_handler` verbatim:

  ```python
  win.tk.eval(f'Ng_LoadMesh "{path}"')   # native loader
  win.tk.eval('set selectvisual mesh')   # *** required ***
  win.tk.eval('Ng_SetVisParameters')
  win.tk.eval('redraw')
  win.tk.eval('Ng_ReadStatus')
  ```
  inside `win.after(100, _load)` and blocks the main thread on
  `win.mainloop()`.  Two earlier wrapper drafts failed:
  `ngsolve.Mesh()`+`Draw()` (NGSolve's Tk-side draw is webgui-only,
  per PR #241 author note) and bare `Ng_LoadMesh` without
  `selectvisual mesh` (mesh loads into global scene but viewer panel
  stays in geometry mode).

* `.sol` support: NGSolve `gfu.Save()` writes a raw float64
  coefficient dump, NOT Netgen's text "solution" format -- so
  `Ng_ImportSolution` cannot read it.  The viewer instead reconstructs
  the FES by matching ndof against `sol_size/8` (real) or `/16`
  (complex), in order: HDiv o=1/2, HCurl o=1/2, H1 scalar o=1/2,
  H1 dim=3 o=1/2.  `gfu.Load(sol)` + `Draw(gfu, mesh, name)` +
  `set selectvisual solution` overlays the field on the mesh.

* Companion `.vol` discovery for `.sol`: peels trailing `_<word>`
  segments from the stem and tries both `<stem>.vol` and
  `<stem>_fem.vol` at each level (Radia panel-output convention:
  `<base>_<solver>_fem.vol` mesh paired with
  `<base>_<solver>_B.sol` / `_Jsurf.sol` / `_qsurf.sol` / `_Jvec.sol`).

Usage:
```
radia-vol-viewer --register      # install .vol/.sol file associations
radia-vol-viewer model.vol       # CLI invocation (stdout)
radia-vol-viewer field.sol       # auto-detects companion .vol + FES
```

The pre-existing top-level `tools/vol_sol_viewer.py` becomes a thin
backwards-compat shim that delegates to `radia.tools.vol_sol_viewer`.

## 4.63.0 — Thermal Method split into 3 explicit choices

Released 2026-05-21.

The single ``Method = "Thermal"`` choice (v4.59.0-v4.62.0) is
replaced by three explicit Thermal solvers in the radia-ih
Method dropdown:

* ``Thermal: 3D static (no rotation)`` -- 3D heat equation,
  rotation_rpm = 0, q_surf held azimuthally fixed.  Use for
  static one-shot heat-up, feasibility studies, non-rotating IH.
* ``Thermal: 3D + rotation (q_surf re-sampled per step)`` -- 3D
  heat equation, workpiece spins around +z axis, q_surf
  re-projected on the body frame each timestep (v4.58.0+
  feature).  Use for non-axisymmetric workpieces OR
  non-axisymmetric coils under rotation.
* ``Thermal: 2D axisymmetric (rotation implicit)`` -- axisym
  (r, z) solver, 10-100x faster than equivalent 3D.  Use for
  rotationally-symmetric workpieces (cylinder, stepped shaft,
  disk).

The Method dropdown owns the (mesh_type, rotation_rpm) pair;
the embedded HeatPanel's individual ``Mesh type`` combo and
``Rotation [rpm]`` line are HIDDEN because the parent Method
already encodes those.  HeatPanel still exposes the choice
internally for direct subprocess invocation (``calc_heat.py`` /
``calc_heat_axisym.py``).

### Breaking changes

* The constant ``radia.radia_ih.METHOD_THERMAL`` is removed.
  Use one of ``METHOD_THERMAL_3D_STATIC`` /
  ``METHOD_THERMAL_3D_ROTATING`` / ``METHOD_THERMAL_AXISYM``,
  or the ``THERMAL_METHODS`` frozenset for membership checks.
* Old scripts that did ``panel._method_combo.setCurrentText(
  METHOD_THERMAL)`` fail with ``NameError``.  Pick one of the
  three explicit methods.

panel_qa.py registry: ``ih_thermal`` replaced by three
entries: ``ih_thermal_3d_static`` / ``ih_thermal_3d_rotating`` /
``ih_thermal_axisym``.

## 4.62.0 — remove radia-heat standalone (HeatPanel moved to _heat_panel)

Released 2026-05-21.

The pre-4.59.0 standalone ``radia_heat.py`` module and its
``radia-heat`` console-script entry point are now fully removed.
Heat analysis lives exclusively as the Method = "Thermal" choice
in ``radia-ih`` (since v4.59.0).  The HeatPanel sub-widget moved
to ``src/radia/_heat_panel.py`` as an internal implementation
detail of the IH panel.

### Breaking changes

* ``from radia.radia_heat import HeatPanel, HEAT_SRC_SPATIAL, ...``
  fails with ``ModuleNotFoundError`` -- update to
  ``from radia._heat_panel import HeatPanel, HEAT_SRC_SPATIAL, ...``
* The ``radia-heat`` console script is removed.  Old shortcuts
  that launched it must be updated to launch ``radia-ih``
  instead, then pick Method = "Thermal" (and fill qsurf .sol +
  em .vol + wp .vol manually, or use the "Run thermal..." chain
  shortcut after a successful EM run).
* In 4.59.0-4.61.0 the ``radia-heat`` CLI was a deprecation stub
  that redirected to ``radia-ih`` automatically; that transitional
  redirect is gone.

Docs (``docs/IH_THERMAL_WORKFLOW.md``) and the
``radia_mcp.ih.thermal`` topic both reflect the removal.

## 4.61.0 — fix rect_united section-plane regression (task #30 close-out)

Released 2026-05-20.

Three independent bugs in ``_filaments_from_section_planes`` that
the build123d -> OCP shim (v4.56.0) surfaced are all fixed.  The
build123d era happened to mask them via offsetting effects;
post-shim the regression appeared as a ValueError or 32% L error
on the rect_torus_lofted_united test fixture.

* **Bug 1**: ``centroids_attempted`` contract violation.  The caller
  passed ``path_m[i]`` (spine point at bbox-heuristic R) where
  ``_filaments_from_per_station_faces`` requires ``face.center()``;
  ``cad_to_m = span_m / span_cad`` then recovered (R_spine /
  R_face) ~ 0.928 instead of 1.0, silently shrinking UV samples
  by ~8%.  Fix: pass ``face_center_m = face.center() /
  cad_units_per_meter``.
* **Bug 2**: ``is_open`` dispatch routed swept-cross-section coils
  to the rim tracer (right for straight LEADS, wrong for gapped-
  arc geometries with cap_a/cap_b).  Fix: when caps are detected,
  use ``_gen_spine`` (planar long-arc) instead of the rim tracer.
* **Bug 3**: ``_parallel_transport_frame`` returns ``v_hat`` aligned
  with the chord-perpendicular (~9 deg off the true radial
  direction at n_stations=20), shrinking the rect's radial extent
  by cos(9 deg).  Fix: detect planar-arc spines and override
  ``v_hat = radial direction`` at each station (non-planar spines
  still use parallel transport).

Test: ``tests/panels/test_inductance_golden.py::
test_inductance_peec_vacuum_rect_united`` now PASSES (was xfail
2026-05-20).  Golden L updated 145.3 -> 191.5 nH; the pre-fix
145.3 was a GEOMETRIC COINCIDENCE (filaments inadvertently at
bbox-spine R=45.9 instead of conductor R=50, but happened to
land near BEM-A's 153 nH).  Post-fix filaments sit at the loft's
actual interior face centroids -- the corrected-PEEC result is
locked, NOT the physically-most-accurate value (PEEC-perimeter
on rect cross-sections has known gaps vs BEM-A surface RWG
which captures rect-corner current crowding).

Production round-wire pipeline (3turncoil, ih_peec_inductance,
gapped_torus) is COMPLETELY UNAFFECTED: those geometries route
through ``_filaments_from_circle_edges_per_station``, which is
not touched by this fix.  Verified L_coil unchanged on both
round-wire fixtures.

## 4.60.0 — T .sol always saved + IH Summary closes the thermal loop

Released 2026-05-20.

Two changes that close the radia-ih thermal pipeline so users get
to the temperature distribution WITHOUT leaving the panel:

### feat(heat): T .sol always saved for later evaluation

calc_heat.py / calc_heat_axisym.py now ALWAYS save the final
temperature GridFunction as a NGSolve ``.sol`` file, symmetric
with the qsurf.sol contract on the EM side.  Two file-layout
branches:

* ``--msh-output FILE.msh`` set: writes ``FILE_T.sol`` +
  ``FILE_heat.vol`` alongside the GMSH bundle.  (Already the case
  pre-v4.60; this release just exposes the paths in JSON.)
* ``--msh-output`` not set: writes ``<wp-stem>_heat_T.sol`` next
  to the input ``wp.vol`` -- re-use ``wp.vol`` as the companion
  mesh.  (New in v4.60; previously this branch produced no .sol.)

JSON additions in both solvers:

  "T_sol_file":    absolute path to the saved T .sol
  "heat_vol_file": absolute path to the companion .vol (empty when
                   the input wp.vol is the companion)

Reload pattern (same 3-rule contract as qsurf.sol):

```python
mesh = Mesh("workpiece_thermal.vol")     # or <msh-stem>_heat.vol
gfT  = GridFunction(H1(mesh, order=1))   # MUST match solve --fes-order
gfT.Load("workpiece_thermal_heat_T.sol")
T_at = float(gfT(mesh(x, y, z)))         # sample at body point
```

Documented in ``docs/IH_THERMAL_WORKFLOW.md`` (new "Reloading the
T .sol" subsection) and in ``radia_mcp.ih.thermal`` topic (new
"Output: T .sol re-loadable for later evaluation" section).

### feat(radia-ih thermal): IH Summary + Open GMSH auto-fire

radia_ih._on_finished now pretty-prints a "--- Thermal ---" block
in the IH Summary when calc_heat / calc_heat_axisym JSON emits
T_max_C / T_min_C / Q_input_J / probe history / msh_file /
T_sol_file / heat_vol_file / vtu_files / rotation_rpm.  Block
includes Delta T from initial, time span + step count, rotation
status, probe start->end trajectory, and output file paths
(T .sol, heat .vol, GMSH .msh, VTU sequence).

Auto-launch GMSH on the T-distribution .msh after a successful
thermal run -- the same call as a manual "Open GMSH" click but
without the extra user step.  EM methods unchanged (button still
enables on msh_file present but does NOT auto-open).

Combined with the prior v4.58.0 (rotation, .sol strict contract)
and v4.59.0 (Method-dropdown integration), the radia-ih workflow
is now CLOSED for the full IH pipeline:

  1. Pick EM method, click Run -> qsurf.sol + em_vol written.
  2. "Run thermal..." button -> Method switches to Thermal +
     qsurf.sol / em_vol pre-filled.
  3. Set wp.vol + material + dt + t_end + rotation, click Run.
  4. T_max / Delta T appear in IH Summary; GMSH auto-opens on the
     T-distribution .msh; T .sol + companion .vol path also
     shown for downstream re-load.

No external tool invocations beyond the auto-opened GMSH viewer.

## 4.59.0 — heat analysis integrated into radia-ih (radia_heat standalone deprecated)

Released 2026-05-20.

Heat analysis is now a Method choice in the radia-ih panel, alongside
PEEC-Inductance / PEEC-BEM / FEM-Kelvin / FEM-coilmesh.  The
HeatPanel sub-widget (from radia_heat.py) is embedded as a section in
IHPanel and becomes visible only when ``method == "Thermal"`` is
selected.  One window for the full EM->Thermal chain.

### User-visible changes

* New method entry **"Thermal (heat transfer from saved q_surf .sol)"**
  in the radia-ih Method dropdown.  Selecting it hides the EM-side
  sections (Drive, Coil material, Coil geometry, Workpiece material,
  Workpiece impedance, Linear solver, Advanced) and shows the
  Thermal sub-panel (mesh type / heat source / qsurf .sol + em .vol
  inputs / material / convection BC / time scheme / probe /
  rotation_rpm).
* The **"Run thermal..."** chain button (active after a successful
  EM solve produced ``qsurf.sol``) now SWITCHES the method dropdown
  to Thermal and pre-fills the embedded HeatPanel's qsurf_sol /
  em_vol fields.  No new window opens; the heat analysis runs in
  the same IHWindow on the next Run click.

### Deprecation: ``radia_heat`` standalone

``radia_heat.py main()`` is now a deprecation stub.  When invoked it
emits a ``DeprecationWarning`` and redirects to ``radia-ih`` with the
Thermal method pre-selected and the CLI flags applied as panel
pre-fills.  Old shortcuts and shipped CLI integrations keep working
during the deprecation window.  The standalone HeatWindow class
remains importable but is **scheduled for removal in the next minor
release**.

### Architectural notes

* IHPanel.is_runnable, .build_command, .wp_vol_path delegate to the
  embedded HeatPanel when method == Thermal.  HeatPanel's own
  build_command picks calc_heat_axisym.py vs calc_heat.py based on
  the mesh_type combo, enforces the v4.58.0 strict
  ``.sol + .vol`` contract, and inherits the v4.58.0 rotation
  feature (``--rotation-rpm > 0`` re-projects q_surf each timestep).
* Previously-unkeyed sections "Drive", "Coil material", "Linear
  solver", "Advanced" now have ``_sec_drive`` / ``_sec_coil_mat`` /
  ``_sec_solver`` / ``_sec_advanced`` keys so they collapse cleanly
  when method == Thermal.

### Test coverage

* ``tests/panels/test_heat_rotation.py``: 3 unit tests still PASS
  (rotation projection math validated against synthetic unit cube;
  unaffected by the integration since it tests calc_heat directly).
* ``tests/panels/panel_qa.py`` registry: ``heat_3d`` / ``heat_axisym``
  standalone HeatWindow entries replaced with a single ``ih_thermal``
  entry that exercises the integrated Thermal method through
  IHWindow.  Existing panel-rendering / font / layout checks now
  cover the integrated path.
* ``tests/panels/test_heat_chain_golden.py``: e2e chain test
  invokes ``calc_heat`` via subprocess directly -- unchanged
  contract, still PASSes.

### Breaking changes

* The "Run thermal..." button no longer launches a separate process.
  Workflows that programmatically detected the radia_heat window are
  now exercised through IHWindow with method=Thermal.
* The standalone ``radia-heat`` CLI continues to work but logs a
  DeprecationWarning at startup.

## 4.58.0 — 3D thermal workpiece rotation + tightened qsurf/em-vol contract

Released 2026-05-20.

### feat(heat): 3D workpiece rotation

`calc_heat.py` previously treated `--rotation-rpm` as metadata only;
q_surf was held azimuthally static so a spinning workpiece would not
actually see a moving heat source.  This release implements true
rotation for the spatial-qsurf path:

* New helper closure in `_build_qsurf_cf`: at angle theta_rad it
  re-samples q_em on the wp surface with the body rotated by
  theta_rad around the z axis -- world coords
  `(x*c - y*s, x*s + y*c, z)` at body point `(x, y, z)`.  The
  GridFunction backing q_cf is updated **in place**; mesh / FES /
  stiffness / mass remain fixed (LinearForm RHS already reassembles
  each step for the convection term, so per-step overhead = one
  re-projection ~10 ms).
* Time loop calls the resampler at start of each step with
  theta = omega_mech * t when `rotation_active` (spatial qsurf +
  positive rpm).  q_int (integrated heat input) is re-computed when
  rotation is active so Q_input_J accounting stays honest if the
  hotspot partially leaves the workpiece face at some angles.
* `--rotation-rpm > 0` with `--q-uniform` is intentionally a no-op
  (constant in space, nothing to rotate); the panel warns.
* `radia_heat.py` "Rotation [rpm]" tooltip updated to describe the
  new behaviour.

Unit tests in `tests/panels/test_heat_rotation.py` validate the
projection math on a synthetic unit cube with q_em(x,y,z)=x: at
theta=0 / pi / pi/2 the body point near (+1,0,0) reads q ~ 1 / -1 / 0.

### fix(heat): --em-vol required when --qsurf-sol is supplied

Tighten the .sol + .vol contract per "No Fallbacks" policy.  The
previous auto-locate of the sibling `<stem>_fem.vol` is removed.
NGSolve .sol is a raw coefficient vector (no embedded mesh, no
fes_order header), so requiring the companion .vol explicitly is
the only safe contract.

Both `calc_heat.py` and `calc_heat_axisym.py` now raise ValueError
when `--qsurf-sol` is given without `--em-vol`.  Error message
points the user at the typical `<stem>_fem.vol` companion.
`radia_heat.py` panel: "EM .vol (auto):" label dropped (now just
"EM .vol:"), tooltip explains the contract, `is_runnable()`
requires both paths to exist, `build_command()` raises if either
is empty (defence in depth).

`radia_ih.py`'s "Run thermal..." launch button still auto-resolves
the sibling and passes BOTH paths explicitly to `radia_heat.py` --
unaffected.

### Breaking changes

* `calc_heat.py --qsurf-sol` users who previously omitted
  `--em-vol` and relied on the `<stem>_fem.vol` auto-locate must
  now pass `--em-vol` explicitly.

## 4.57.0 — proximity-aware iterative PEEC (default ON in calc_inductance panel)

Released 2026-05-20.

New module `radia.peec_proximity.solve_proximity_iterative` wraps the
existing perimeter PEEC loop-bundle solve with an outer iteration
that augments per-filament `Zs_fil` with the Leontovich surface
dissipation `½ Re(Z_s) · |H_t|² · s_k · L_k` evaluated from the
actual Biot-Savart H field at each filament's wire-surface position.
The bundle re-solves; `Zs_fil` updates with under-relaxation
(default relax=0.3) until `‖ΔZ‖/‖Z‖ < tol` (default 1e-3).

Validated against the 2-wire analytical (50.5 / 23.2 / 12.1 A/m at
three eval positions, exact match).  Critical kernel details
(debugged 2026-05-19/20):

- Biot-Savart angular sign convention: `l·r1/|r1| - l·r2/|r2|`
  (cos α₁ - cos α₂ form).
- r_min clamp must be **1e-12**, NOT 0.1·wire_radius (the latter
  over-clamps close-range contributions, giving |H_t_self| ~360×
  too small).
- Self-exclusion is per-segment only; adjacent segments of the same
  filament must contribute to |H_t| at the eval surface point.

Production result on the 3-turn pancake (Cu, 150 kHz, n_peri=16):

| Path | R | L |
|---|---|---|
| pre-2026-05-19 (R_DC only) | 0.3945 mΩ | 426.30 nH |
| 2026-05-19 v4.55.4 (Bessel self-skin) | 3.6752 mΩ | 430.14 nH |
| **v4.57.0 (+ proximity)** | **4.4793 mΩ** | **431.32 nH** |

`proximity_factor = 4.4793 / 3.6763 = 1.218` at 150 kHz; the factor
is constant ~1.20-1.22 across 1 kHz - 1 MHz on this geometry — the
**structural ceiling of perimeter PEEC**, which captures
surface-Leontovich proximity but **not** transverse eddy loops in
the wire interior.  The remaining gap to the 15 mΩ LCR hi-tester
measurement is most plausibly lead/contact resistance (5-50 mΩ
typical for 2-terminal probes at low Ω); 4-terminal Kelvin
re-measurement is the deciding test.  Volume PEEC (radial filaments
beyond perimeter) is sketched in `docs/peec/VOLUME_PEEC_DESIGN.md`
and marked DEFERRED pending that re-measurement.

`calc_inductance.py` exposes the new `--peec-proximity` /
`--no-peec-proximity` flag (default ON).  The opt-out restores
v4.55.4 self-only Bessel behaviour (R=3.6763 mΩ on the 3-turn case
to within numerical noise).

Golden tests updated:
- `tests/panels/golden/peec_inductance_3turn_150kHz_Cu.json`:
  L 426.30 -> 431.32 nH, R 0.3945 -> 4.4793 mOhm (both rate-limiters
  through the two-step lift R_DC -> Bessel -> Bessel+proximity).
- `tests/panels/golden/peec_inductance_torus_50kHz_Cu.json`:
  L 85.10 -> 87.08 nH, R 0.118 -> 0.6496 mOhm (this golden was last
  touched in v4.12.0 and missed the v4.55.4 Bessel fix; sync'd here).

Bundled performance fixes (Option B prerequisites):
- `_b3d_shim.Solid`: lazy caches for bounding_box / faces / edges /
  solids.  `BRepBndLib.AddOptimal_s` is ~0.5 s on the 3turncoil
  sample and was being re-computed 6× by independent validation
  predicates in `coil_from_cad.filaments_from_step`.
- `coil_from_cad._filaments_from_circle_edges_per_station`:
  vectorise the O(N·K) centre-proximity dedup with
  `numpy.norm(axis=1).min()`.  Drops Python any()/norm() count
  from ~199k to ~5k on the same sample (~0.9 s recovered on the
  9.4 s cold run).

Pre-existing unrelated bug (filed as separate task): `_filaments_from_section_planes`
fails on `tests/panels/golden/rect_torus_lofted_united.step` with
"spacing between cap_a and the adjacent interior station is
21.96 mm vs typical 4.60 mm".  Introduced when the OCP shim landed
in v4.56.0; the PEEC `rect_united` golden has been failing on main
since then.  Out of scope for v4.57.0; will be fixed separately.

## 4.55.0 — endpoint anchoring to cap centroids (eliminates rim-end kink at lead caps)

Released 2026-05-16.  Final fix in the keiko viz response chain.
v4.54.0 reduced corner bunching via Wang-Joe RMF + densification but
left a residual 41.6 deg kink at vertex N-1 visible in GMSH near the
cap tips (keiko's 2nd viz report image).  Root cause: the "longest
open edge" spine traces the conductor's LATERAL RIM (z=+wire_radius
on a flat coil), while interior centerline points correctly come from
midpoint sectioning -> face centroids (z=0).  The path_cad endpoints
path_cad[0] / path_cad[-1] were pinned to the rim endpoints, creating
a kink at the rim-to-centroid transition in the final segment.

Symptoms (keiko's outsideline.step):
- 41.6 deg bend at vertex N-1 (visible kink in GMSH near each cap)
- 48% spread in |I| across the 16 perimeter filaments, distorted by
  the kink coupling asymmetrically into L

Fix: replace the rim endpoints with the cap-face centroids from
`coil_topology.extract_coil_topology(solid).cap_a/b.center()`.  Map
each spine endpoint to its NEAREST cap centroid (handles the cap_a /
cap_b ordering automatically).

Results after fix:
- Max bend per step: 41.6 deg -> 24.9 deg (cap-end kink eliminated;
  the residual 24.9 deg is the TRUE structural lead-arc 64 deg corner
  densified by v4.54.0 to 24.9 deg per step)
- Zero bends > 30 deg (was 1)
- L_coil 91.82 -> 92.22 nH (+0.4%, the corrected geometry slightly
  changes the L matrix; the new value is more physically accurate)
- |I| distribution now reflects the TRUE inner-vs-outer asymmetry
  of a curved wire (wire_r=3mm, arc_R=30mm, ratio 0.1 -> ~55% spread
  expected from physics).  Previous distorted distribution was an
  artifact of the rim-end kink.

### Coordinated bumps

- radia 4.54.0 -> 4.55.0
- radia-mcp 0.54.0 -> 0.55.0
- cubit-mesh-export unchanged

55 passed, 1 skipped on the coil-pipeline regression suite (same
count: endpoint anchoring is a behavior change, not a new test).

## 4.54.0 — RMF (Wang-Joe) + corner densification for filament viz smoothing

Released 2026-05-16.  Follow-up to v4.53.0 addressing keiko's viz
report (filament "bunching" visible at lead-arc junction in GMSH).

### Background

Even though v4.53.0's adaptive resampling gives correct PEEC L_coil
(keiko's STEP -> 91.70 nH, physical), GMSH visualization showed the
16 perimeter filaments compressing at the corner.  Diagnosis: this
is partly FORESHORTENING (the cross-section plane rotates with the
spine bend; viewed perpendicular to the bend axis it projects to a
narrower visual extent) and partly POLYLINE COARSENESS (each
filament's polyline has a single hard kink at the corner, no
smoothing).

### What shipped

**Rotation-Minimizing Frame (RMF) replaces parallel-transport**:
`_parallel_transport_frame` now implements Wang-Joe double-reflection
(Wang et al. 2008, "Computation of rotation minimizing frames",
ACM TOG 27(1):2) instead of incremental Rodrigues rotations.  The
double-reflection method minimizes the integral of squared angular
velocity along the curve -- provably-minimum accumulated twist on
polylines with kinks.  Function name retained for backward compat.

On the keiko fixture, RMF gives L_coil = 91.82 nH vs v4.53.0
parallel-transport's 91.70 nH -- 0.13% difference, within numerical
noise (both methods agree on circular cross-section; only differ on
twist accumulation which is a secondary effect).  Tests added:
- `test_rmf_orthogonality`: strict (t, u, v) orthonormal at every
  vertex including sharp kinks.
- `test_rmf_twist_minimization_vs_pt_on_straight_path`: zero twist
  on straight path (frame must NOT rotate).

**Corner densification** in `_centerline_from_open_spine`: after
adaptive resampling, insert intermediate spine points ON THE OCC
CURVE near sharp bends until the per-segment bend angle <=
`max_bend_per_step_deg` (default 20 deg).  Each insertion bisects
the LARGER of the worst vertex's two neighbours to balance segment
lengths.  Hard cap at 500 total points to prevent infinite blowup
on numerical jitter.

On the keiko fixture, the spine densifies from 68 -> 80 stations
(+12 inserted near the lead-arc corners), reducing max bend per
step from 64 -> 41.6 deg with mean bend 6.3 deg / median 5.6 deg.
The 41.6 deg residual is the TRUE polyline kink at the lead-arc
junction -- bisecting a linear segment around it gives the same
bend at the bisection midpoint, so the kink itself cannot be
softened.  But the SURROUNDING bends are smoothed, which is what
the viz needs.

### Performance

RMF: ~same cost as PT (both O(N) per vertex).  Corner densification:
adds a few iterations of segment scanning + OCC `spine @ t` queries
near corners -- typically +5-20 stations on coil geometries, ~10-50
extra ms.

### Coordinated bumps

- radia 4.53.0 -> 4.54.0
- radia-mcp 0.53.0 -> 0.54.0
- cubit-mesh-export unchanged

55 passed, 1 skipped on the coil-pipeline regression suite (was 52/1
in v4.53.0: +3 RMF + densification tests).

## 4.53.0 — keiko's "arc + leads" 1-turn coil now works end-to-end (CCW winding + adaptive resampling)

Released 2026-05-16.  Responds to keiko's patch report
(`W:\31_Go-Tech\10_IH_toymodel\2026_05_10_1turn_coil_loft\README_coil_from_cad_patch.txt`)
which documents 3 patches she applied locally to get her
`1turn_coil_loft_outsideline.step` PEEC-solvable.  v4.53.0
integrates her bug fix verbatim (CCW winding) and replaces her
spine-thinning workaround with a design-level adaptive sampling
that is policy-compliant (CLAUDE.md "No Fallbacks").

### Fix 1: CCW winding in `_sample_face_perimeter_in_pt_frame`

**Bug**: per-segment Cubit lofts produce shared cross-section faces
whose outer-wire orientation alternates between adjacent volumes.
The unmodified sampler returns CW samples at alternating stations;
the downstream parallel-bundle solver assumes consistent orientation
and connects sample k at station i to sample k at station i+1.
With orientation flip, adjacent stations zigzag in opposite
directions around the cross-section -- filament paths self-intersect
and the Ruehli L matrix degenerates to NaN.

**Fix** (keiko's patch verbatim): after projecting samples to (u, v)
in the parallel-transport frame, compute the signed polygon area
(`sum(u[i]*v[i+1] - u[i+1]*v[i])`) and reverse the array if
negative.  Idempotent on already-CCW samples.

### Fix 2: Adaptive resampling in `_centerline_from_open_spine`

**Background**: keiko reported that hardcoded `n_segments=100` on
her 208 mm spine gives 2.08 mm segments -- with 2.9 mm wire radius
the ratio is 0.72 (below the 1.0 threshold), and the
lead-arc 64 deg bend tripped `_check_spine_no_singular_corner`
even though the geometry is physically fine.

**Her proposed workaround** (Fix 2/3 in the README): iteratively
delete the shortest segment until ratio >= 1.  This is a SILENT
geometry modification -- conflicts with CLAUDE.md "No Fallbacks"
(she acknowledged this in the README).

**v4.53.0 alternative**: adaptive resampling UPFRONT.  Before the
spine sampling loop, probe one OCC section at the midpoint to
estimate the wire radius, then cap `n_segments` at
`floor(spine_length / (1.10 * wire_radius_estimate))`.  Result:
keiko's 208 mm spine resamples to 65 segments (3.2 mm each, ratio
1.10), which clears the singular-corner check by design.  The
check still fires on genuinely degenerate geometries (caller-
supplied path with literal overlap, etc.).

This is NOT a fallback -- the design contract is "adjacent stations
are at least 1 wire-radius apart so parallel transport is stable".
The caller's `n_segments` is interpreted as an upper bound, not an
absolute count (the function already takes liberties with the
caller's request to satisfy other constraints).

### Test: keiko's STEP now produces L_coil = 91.70 nH

Replaced the v4.49.0 `test_keiko_outsideline_filaments_raises_on_singular_L`
with `test_keiko_outsideline_succeeds_with_adaptive_resampling`,
asserting `np.isfinite(L).all()` and `50 < L_coil_nH < 300`.

On the synthetic vertex-aligned replica (Predicate 1 UV path,
unchanged from v4.49.0), L_coil = 142 nH.  On keiko's actual
non-vertex-aligned STEP (Predicate 4 + adaptive resampling),
L_coil = 91.70 nH.  Different by 35% because the spine traces a
slightly different path through the lead-arc junction; both are in
the physical range for a 1-turn 30 mm-radius coil + leads.

### Coordinated bumps

- radia 4.52.0 -> 4.53.0
- radia-mcp 0.52.0 -> 0.53.0
- cubit-mesh-export unchanged

52 passed, 1 skipped on the coil-pipeline regression suite (same
count as v4.52.0; -1 raises-on-singular test, +1 succeeds-with-
adaptive test).

## 4.52.0 — final magic-number pins + negative-confidence tests

Released 2026-05-16.  Closes the PEEC STEP-loading magic-number
audit by pinning the remaining load-bearing factors that the
v4.49.0/v4.50.1 boundary tests did not cover.  No functional
changes; +6 tests.

### Magic-number pins added (4)

- `test_dedup_tol_circle_edges_pinned_at_0_1_median_r`
  (`coil_from_cad.py:843` / `:1064`, semicircle pair dedup tolerance
  on Cubit unite-split CIRCLE edges)
- `test_dedup_tol_loft_cross_sections_pinned_at_0_1_eq_radius`
  (`coil_from_cad.py:1205`, shared-end-cap merge on NON-united loft
  STEPs)
- `test_detect_lead_bars_radius_spread_pinned_at_0_1`
  (`coil_from_cad.py:497`, CYLINDER face radius / median tolerance
  for lead-bar acceptance)
- `test_detect_lead_bars_length_factor_pinned_at_5_0`
  (`coil_from_cad.py:515`, lead length >= 5.0 * wire_radius
  threshold)

### Negative-confidence tests added (2)

- `test_predicate_3_does_not_fire_on_keiko_no_revolution_faces`:
  keiko's 2-BSPLINE + 2-PLANE STEP must NOT match Predicate 3
  (which requires a revolution-type face); dispatch must reach
  Predicate 4.
- `test_predicate_5_does_not_fire_on_open_coil`: keiko's STEP has
  2 cap faces so `extract_coil_topology` returns `is_open=True`;
  Predicate 5 (CLOSED-only) must NOT fire.

### Coordinated bumps

- radia 4.51.0 -> 4.52.0
- radia-mcp 0.51.0 -> 0.52.0
- cubit-mesh-export unchanged

52 passed, 1 skipped on the coil-pipeline regression suite (was
46/1 before this release: +6 pin + negative tests).

### PEEC STEP-loading sweep -- audit complete

Across v4.48.2 / v4.49.0 / v4.50.0 / v4.50.1 / v4.51.0 / v4.52.0,
every layer of the pipeline now has fail-fast positive proof:

| Layer | Check |
|-------|-------|
| CAD entry | multi-solid raise |
| Centerline dispatch | classification single dispatch |
| Centerline result | bbox-containment + per-point distance |
| Filament construction | corner detect (pre-Ruehli, covers HACApK) |
| Filament topology | bbox-cover + inside-bbox + near-surface |
| Solver assembly | finite-L safety net (post-Ruehli) |

19 load-bearing magic numbers pinned by boundary/source tests.

## 4.51.0 — Strong Tier C: per-point distance-to-solid via BRepExtrema

Released 2026-05-16.  Strengthens the v4.50.0 centerline-inside-bbox
check with an additional STRONG positive proof: each centerline
point's distance to the solid boundary (via `BRepExtrema_DistShapeShape`)
must not exceed `1.10 * wire_radius`.  This catches the failure modes
that bbox-containment misses:

- Predicate 4 picking a surface-rim edge whose centerline lies ON the
  solid surface but INSIDE the bbox.
- Predicate 5 racetrack-as-circle where the circle's diameter is
  inside the racetrack bbox but the circle exits the conductor's
  cross-section at the corners.
- Wrong-radius spine that stays within bbox but escapes the wire tube.

### Design notes

- **BRepExtrema_DistShapeShape works on BSpline solids** (verified
  2026-05-16 against the v4.50.0 `BRepClass3d_SolidClassifier`
  reliability issue): points INSIDE return d=0, points OUTSIDE return
  the actual distance to the nearest surface.  This is a per-point
  computation, not a topological classification, so the BSpline
  classifier's geometric edge-case failures do NOT apply.
- **Subsampling**: BRepExtrema is O(face_count) per point, ~17 ms/pt
  on a 4-face STEP -> ~17 s on 100 pts of a 700-face STEP.
  Sub-sampling to 20 evenly-distributed points keeps the check
  bounded at ~2 s on a 700-face STEP.  Sub-sampling is sufficient
  because the wrong-location failure modes affect contiguous regions,
  not isolated points.
- **Distance tolerance 1.10 * wire_radius**: empirically 100% of
  per-station-mean centerlines on a smooth build123d sweep coil fall
  within `wire_radius` of the lateral surface (parallel-transport
  displacement is bounded).  10% slack covers numerical noise on
  the tube boundary.
- **Three-layer orthogonal verification** (NOT a fallback chain):
  bbox-cover (filaments under-coverage) + inside-bbox (gross
  wrong-location) + near-surface (per-point envelope).  All three
  must pass.

### Tests added (4)

- `test_check_near_solid_surface_accepts_wire_axis`: cylinder axis
  centerline passes.
- `test_check_near_solid_surface_rejects_far_off_axis_spine`:
  spine 10mm outside a 5mm-radius cylinder fails (= 10mm exit > 5.5
  mm tolerance).
- `test_check_near_solid_surface_distance_tolerance_pinned`: pins
  the 1.10 factor against drift.
- `test_predicate_1_does_not_fire_on_keiko_split_lateral`: negative
  confidence -- the keiko 50/50 split-lateral STEP must NOT match
  Predicate 1, so dispatch reaches Predicate 4 where the
  singular-corner check fires.

### Coordinated bumps

- radia 4.50.1 -> 4.51.0
- radia-mcp 0.50.1 -> 0.51.0
- cubit-mesh-export unchanged

46 passed, 1 skipped on the coil-pipeline regression suite (was 42/1
before this release: +4 strong-Tier-C + negative-confidence tests).

## 4.50.1 — doc lies cleanup + readability + magic-number pin tests

Released 2026-05-16.  Polish pass after the v4.48.2/v4.49.0/v4.50.0
PEEC STEP-loading sweep.  No functional changes; addresses the
remaining review smells.

### Doc lies removed (cascade-era residue)

- `filaments_from_step` docstring no longer claims "Falls back to
  legacy path if CoilBuilder reconstruction fails" -- the two paths
  are SEPARATE caller-chosen entries per CLAUDE.md "No Fallbacks";
  failures propagate.
- `_centerline_from_open_spine` docstring no longer says "Caller
  should fall back to `_centerline_from_cross_sections`" -- a failure
  here means the dispatcher's positive predicate was wrong, not that
  another extractor should be tried.
- `_collect_circle_edge_centers` docstring no longer claims "lets
  the caller fall through to other extractors" -- returning None is
  a predicate-style negative match, not a fallback signal.
- `_check_filaments_cover_solid_bbox` setup comment rewritten to
  describe the orthogonal-pair check policy (bbox-cover + inside-
  bbox), removing the historical "v4.38.0 Tier-2c check did NOT
  cover Tier 1/2/2b" wording.

### Readability

- `peec_bundle.py:176` `diag_bad` calculation replaced the
  convoluted `int(np.isfinite(np.diag(L)).sum() - L.shape[0]) * -1`
  with the equivalent `int((~np.isfinite(np.diag(L))).sum())`.

### Magic-number pin tests added (4)

- `test_detect_cap_faces_area_ratio_threshold_accepts_clear_caps`
  (`coil_topology.py:107`, area_ratio_threshold=2.0)
- `test_R_spine_0_85_factor_at_coil_topology_py_149`
  (R_spine = 0.85 * R_outer; load-bearing for Predicate 5 spine)
- `test_check_filaments_cover_solid_bbox_slack_factor_1_5_pinned`
  (slack = 1.5 * wire_radius; load-bearing for under-coverage catch)
- `test_check_centerline_inside_solid_slack_0_05_pinned`
  (slack = 0.05 * bbox diagonal; load-bearing for racetrack catch)

42 passed, 1 skipped on the coil-pipeline regression suite (was 38/1
before this release: +4 pin tests).

### Coordinated bumps

- radia 4.50.0 -> 4.50.1
- radia-mcp 0.50.0 -> 0.50.1
- cubit-mesh-export unchanged

## 4.50.0 — Tier C: universal centerline-inside-bbox positive proof on all 5 predicates

Released 2026-05-16.  Completes the PEEC STEP-loading weakness
campaign started in v4.48.2.  The bbox-cover check
(`_check_filaments_cover_solid_bbox`) has been ORTHOGONALLY paired
with a new `_check_centerline_inside_solid` check that runs on the
centerline (not the filaments) at every extractor's return.  Both
checks must pass -- they catch disjoint failure modes:

- **bbox-cover** (filaments): under-coverage (lead missed, spine
  too short).  Existing since v4.47.x.
- **inside-bbox** (centerline, NEW v4.50.0): wrong-location spine
  (e.g. Predicate 5 mapping a racetrack to a planar circle whose
  diagonal corners exit the racetrack's narrower-axis bbox).

### What shipped

- **`_check_centerline_inside_solid`** in `coil_from_cad.py`:
  verifies every centerline point lies within solid bbox + 5% of
  bbox diagonal slack.  Raises ValueError with the worst-offending
  point coordinates, the bbox extents, and a HINT on supported
  topologies.
- **`_centerline_from_filament_paths`** helper: derives an effective
  centerline from `filament_paths` (per-station mean of n_peri
  filaments) for the Path 1/2/2b/2c entry points that build
  filaments without an explicit centerline.
- **Universal wiring** in `extract_centerline_from_step`: a small
  `_dispatch_and_verify` wrapper runs the inside-bbox check after
  every Predicate 1-5 extractor.
- **Universal wiring** in `filaments_from_step`: a small
  `_verify_topo` helper runs BOTH checks (bbox-cover + inside-bbox)
  after every Path 1/2/2b/2c.

### Design decision: bbox-containment, not SolidClassifier

The original review recommended `BRepClass3d_SolidClassifier` for
a strong "point is inside the BREP solid" test.  Empirical probing
2026-05-16 showed the classifier produces unacceptable
false-positives on BSpline-lateral solids: 78% of true-interior
centerline points on a smooth vertex-aligned sweep coil classified
as OUT (including the wire axis along a straight lead).  This
would have broken every existing BSpline PEEC test case.

The weaker bbox-containment check has:
- ZERO false-positives across the regression suite (38 passed,
  was 35).
- Sufficient strength to catch the racetrack-as-circle Predicate 5
  failure that motivated the review (proven by a synthetic
  boundary test).
- O(N) cost vs O(N * face_count) for the classifier.

A stronger per-point distance-to-surface check is deferred until
a reliable OCC inside-test API surfaces (radia issue tracker
entry pending).

### Boundary tests added

- `test_check_centerline_inside_solid_accepts_centerline_inside_bbox`:
  pins the happy path.
- `test_check_centerline_inside_solid_rejects_racetrack_as_circle`:
  builds a thin rectangular ring + a 28 mm circle centerline,
  asserts ValueError "extends beyond solid bbox".
- `test_check_centerline_inside_solid_slack_accommodates_cap`:
  pins the 5% slack -- cylinder caps with centerline endpoints
  on the bbox extreme must NOT raise.

### Coordinated bumps

- radia 4.49.0 -> 4.50.0
- radia-mcp 0.49.0 -> 0.50.0
- cubit-mesh-export unchanged

38 passed, 1 skipped on the coil-pipeline regression suite (was 35/1
before this release: +3 Tier C boundary tests).

This **completes the PEEC STEP-loading weakness sweep**: Tier A
(silent fallbacks), Tier B (corner detect at construction), Tier C
(centerline positive proof), Tier D (entry guards), Tier E (magic
number tests).  The pipeline is now fail-fast at every layer with
actionable diagnostics, covering both dense and HACApK paths.

## 4.49.0 — Tier A+B+D+E weakness sweep (no-fallback hardening)

Released 2026-05-16.  Multi-front improvement pass after the PEEC
STEP-loading review identified silent fallbacks, layer-violations,
and absent entry guards that the v4.48.x classification dispatch
left in place.  All fixes obey CLAUDE.md "No Fallbacks - Fail
Fast, Fail Loud" -- there is no automatic recovery; the user gets
a hard error with a HINT.

### Tier A: silent fallbacks removed in `_centerline_from_open_spine`

- **`coil_from_cad.py:1255`**: `max(open_edges, key=length)` was
  non-deterministic on equal-length BSPLINE rim halves (Cubit/OCC
  listing-order tie).  Replaced with deterministic
  `(length, centroid_x, centroid_y, centroid_z)` lex sort.
- **`:1278-1281`** `try: section() except: cross = None`: removed,
  raises ValueError with station index, midpoint, tangent, and the
  underlying OCC exception name.
- **`:1294-1297`** "use previous width on section fail": removed.
  Section failure is hard (the open-spine assumption is broken on
  this geometry).
- **`:1295` bug**: `widths_cad[max(0, i-1)]` at i=0 copied 0.0
  (uninitialized).  Now unreachable -- the only path that hit it was
  the silent fallback above.

### Tier B: corner detection moved to filament construction layer

- **New `_check_spine_no_singular_corner`** in `coil_from_cad.py`:
  scans interior spine vertices, raises if `bend > 60 deg AND
  adj_min_segment_length < wire_radius` (the keiko condition,
  captured exactly).  Runs as the first step inside
  `filaments_from_polyline` BEFORE the O(N^2) Ruehli build.
- **HACApK path coverage**: because the check runs before solver
  assembly, both dense-L and HACApK paths are protected by the
  same single check.  v4.48.2's
  `peec_bundle._assert_solver_L_finite` (post-assembly safety net)
  remains as belt-and-suspenders.
- **Diagnostic improvement**: error names the offending spine vertex
  in mm coordinates, the bend angle, and the seg/radius ratio --
  far more actionable than v4.48.2's "non-finite entries at L[395,396]".

### Tier D: entry guards

- **Multi-solid STEP raise** in `extract_centerline_from_step`: a
  STEP containing more than one solid raises ValueError with bbox
  summary (single-coil PEEC handles one solid).
- **`cad_to_m` silent 1.0/0 bug fix** at
  `coil_from_cad.py:_filaments_from_per_station_faces` (line ~2692):
  the previous `norm(c0_m) / norm(c0_cad)` form silently degenerated
  to 1.0 when c0 happened to lie near the origin (quarter-symmetry
  coil with cap at (0, +y, 0)).  Replaced with station-spacing-based
  recovery (`norm(c1 - c0)`), which is origin-independent.  Raises
  if the spacing itself is degenerate (single-station path).

### Tier E: magic number boundary tests

- **New `tests/coil_from_cad/test_predicate_boundaries.py`**:
  - `dominance >= 0.8` accept (clean torus) / reject (keiko 50/50 split)
  - `_check_spine_no_singular_corner` accepts smooth (30 deg, long
    segs), accepts below-threshold (59.9 deg short segs), accepts
    long-seg sharp (90 deg, 5x wire_r), rejects keiko-class
    (64 deg, 0.17x wire_r), rejects zero-length segments
  - multi-solid STEP raises

### What this DOES NOT include (deferred to v4.50.0)

- **Tier C** (universal positive proof, "centerline-inside-solid"):
  requires `BRepClass3d_SolidClassifier` integration + broad-fixture
  validation (rect_torus_lofted_united, ih_closed_torus_coil,
  3turncoil united + non-united, etc.) to confirm no behavior shift
  on existing passing geometries.  Rolling separately so it can be
  reverted if a fixture regresses.

### Coordinated bumps

- radia 4.48.2 -> 4.49.0
- radia-mcp 0.48.5 -> 0.49.0
- cubit-mesh-export unchanged

35 passed, 1 skipped on the coil-pipeline regression suite (was 27/1
before this release: +8 boundary tests).

## 4.48.2 — PEEC bundle solver fail-fast on non-finite L (silent NaN -> ValueError)

Released 2026-05-16.  Hot-fix for the symptom that surfaced AFTER
v4.48.1 routed keiko's `1turn_coil_loft_outsideline.step` through the
correct OPEN longest-edge predicate: the build succeeded, the bbox
covered the leads, but `compute_port_impedance` returned
`L_coil_nH = NaN` -- a "successful" run with garbage numbers (CLAUDE.md
"No Fallbacks - Fail Fast, Fail Loud" violation).

### Why

The OPEN longest-edge path samples the 16 cross-section perimeter
filaments via parallel transport of the spine's tangent frame.  At
the lead-cap junction (a 64 deg corner in keiko's spine), parallel
transport places adjacent segments on the SAME filament near-coincident
(verified: `fil1.seg195 <-> fil1.seg196`, `fil5.seg0 <-> fil5.seg4`,
`fil9.seg195 <-> fil9.seg196`, `fil13.seg0 <-> fil13.seg4`).  The
Ruehli mutual-inductance kernel is singular on coincident segment
pairs and the dense kernel returned NaN / Inf silently.  Downstream
`np.linalg.solve` propagated the NaN into the port impedance.

The user saw `exit_code: 0`, `n_filaments: 16`, `L_coil_nH: NaN` --
the worst possible failure mode (silently wrong "OK" output).

### What shipped

1. **`peec_bundle.build_bundle_solver` now calls
   `_assert_solver_L_finite`** immediately after constructing
   `PEECCircuitSolver`.  If `solver.L` has any non-finite entries,
   raises `ValueError` with: NaN count, diagonal-NaN count, sample
   degenerate filament/segment pairs (first 8), and a HINT pointing
   at the v4.48.1 vertex-aligned-loft fix path (Predicate 1 UV
   sampling, which does NOT use parallel transport).
2. **Regression test rewrite** in
   `tests/coil_from_cad/test_keiko_outsideline_centerline.py`:
   - `test_keiko_outsideline_filaments_raises_on_singular_L`:
     asserts that the keiko fixture now raises ValueError instead of
     returning a NaN-laden topo.
   - `test_vertex_aligned_replica_routes_predicate_1_uv` (new):
     builds a vertex-aligned synthetic STEP via build123d
     `sweep(circle, smooth_spline)`, asserts Predicate 1 UV path is
     hit, asserts `np.isfinite(L).all()`, asserts L_coil in
     [50, 300] nH (physical band for 1-turn 30 mm coil + leads).
3. **`radia-mcp` 0.48.5** coordinated bump.

`cubit-mesh-export` unchanged.

## 4.48.1 — STEP centerline classification dispatch + path_points_m removal

Released 2026-05-15.  Hot-fix for keiko's `1turn_coil_loft_outsideline.step`
on mdx (4 PEEC runs FAIL with "Filament path does not match the
conductor solid bbox").

### Why

The STEP file is a BSPLINE-lofted "arc + 2 lead-bars" coil: 4 faces
(2 BSPLINE lateral halves + 2 PLANE caps), bbox y=[-33, +50] mm with
the leads extending out to y=+50 mm.  The pre-v4.48.1 dispatcher in
`extract_centerline_from_step` was a try/except cascade
`revolution_sweep -> topology_spine -> open_spine`.  The middle path
`_centerline_from_topology_spine` "succeeded" with a planar arc at
R = 0.85 * bbox_max = 42.5 mm that overshot the conductor in x
(12 mm) AND missed the lead bar in y (6 mm), and the longest-edge
fallback `_centerline_from_open_spine` (which produces the CORRECT
spine for this geometry, x=+-30 mm + lead at y=+50 mm) was shadowed.
The downstream `_check_filaments_cover_solid_bbox` sanity check fired
and told the user to pass `--path-points-m` -- but the auto-detect
path that handled their geometry already existed; it was the eager
upstream branch that was wrong.

### What shipped

1. **Classification-based dispatch** in `extract_centerline_from_step`
   (No-Fallback policy, CLAUDE.md "No Fallbacks - Fail Fast, Fail
   Loud").  Five positive-match predicates instead of a try/except
   cascade: multi-station loft -> `_cross_sections`; united multi-
   turn -> `_circle_edge_centers`; revolution+plane -> `_revolution_
   sweep`; OPEN -> `_open_spine`; CLOSED -> `_topology_spine`.
2. **`_centerline_from_topology_spine` is now CLOSED-only** and raises
   if called on an OPEN coil (programming-error indicator, not a
   soft-fallback signal).
3. **`path_points_m` parameter removed** from `filaments_from_step`.
   STEP files are now the single source of truth for the centerline.
   If auto-detection cannot recover a covering spine the user fixes
   the CAD rather than papering over breakage with a hand-crafted
   JSON.  The HINT in `_check_filaments_cover_solid_bbox` was updated
   to drop the `--path-points-m` reference.
4. **Regression fixture + tests** at
   `tests/coil_from_cad/fixtures/keiko_outsideline.step` exercising
   lead-bar coverage, x-bbox no-overshoot, equiv wire radius
   recovery, and the API contract (`path_points_m` removed).
5. **`radia-mcp` 0.48.3** coordinated patch bump per release-triple
   policy.

`cubit-mesh-export` is unchanged this release.

## 4.28.1 — radia_ih Run button stays disabled after Browse... fix

Released 2026-05-08.  Hot-fix for kubota's report on mdx + 100号機.

### Why

When `radia_ih` was launched directly (`python -m radia.radia_ih` /
console entry, no `--vol` arg) and the user picked a `.vol` via the
"Browse..." dialog, the Run button stayed grayed out for every non-
vacuum method (PEEC+BEM, BEM-A+BEM, PEEC+FEM+Kelvin, Full FEM A-V).
Two parallel bugs:

1. `AnalysisWindow._browse_vol()` updated the `.vol` line edit text but
   never re-inspected the file.  `IHPanel._vol_mats` stayed `None`,
   `is_runnable()` returned `False`, the Run button never enabled.
2. `IHWindow.__init__` called `_reload_vol_info(vol_path)` using the
   constructor argument, so even when `_restore_settings()` repopulated
   the line edit from saved `radia_ih.json`, the QSettings-restored
   path was never inspected.

### What changed

* `radia_gui_base.AnalysisWindow` — added overridable
  `_on_vol_changed(path)` hook, invoked from `_browse_vol()` and from
  `QLineEdit.editingFinished` (manual edits).  Default no-op so panels
  that do not depend on `.vol` contents (radia_em, radia_pcb,
  radia_heat) are unaffected.
* `radia_ih.IHWindow._on_vol_changed` — re-runs `_reload_vol_info` +
  `_update_run_state`.
* `radia_ih.IHWindow.__init__` — now passes `self._vol_edit.text()`
  (post-restore) to the initial `_reload_vol_info` instead of the
  stale arg.

Headless QApplication test confirms: fresh launch with empty `vol_path`
leaves Run disabled (expected); Browse... to a valid `.vol` with `sibc`
boundary now enables Run for PEEC+BEM and PEEC+FEM+Kelvin.

## 4.26.0 — BEM-A coil migrated from intree (Python) to ngsolve.bem; intree code retired

Released 2026-05-03.  Strategic pivot after benchmarking.

### Why

The intree BEM-A assembler shipped in 4.25.0
(`radia.bem.efie_rwg.solve_inductance_source_sink_intree`) was a pure-
Python double loop over N² triangle pairs with no C++ acceleration.
Benchmarking on `tests/panels/golden/rect_torus_lofted_united.step`
across N=302..1014 triangles showed ngsolve.bem (`use_fmm=False` +
`mat.COO()` extract) is **50-60x faster at every measured N** with no
crossover.  L values agree to 0.025-0.05 %, so the methods are
numerically equivalent — only speed differs.

| maxh (mm) | n_tris | intree (Python) | ngsolve.bem | speedup |
|-----------|-------|-----|-----|-----|
| 12        | 302   | 23.5 s | 0.39 s | 60x |
| 8         | 340   | 28.8 s | 0.47 s | 61x |
| 6         | 436   | 44.2 s | 0.81 s | 55x |
| 4         | 1014  | 209 s  | 5.91 s | 35x |

### What changed

* `src/radia/bem/coil_inductance_ngsolve.py` — new (promoted from
  `examples/induction_heating/bem_reference/bem_inductance.py`).
  Adds:
    * `compute_inductance_source_sink(mesh, ..., Z_s_re=...)` — returns
      `L`, `R` (AC SIBC), `J`, `gf_J`, etc.  ngsolve.bem LaplaceSL on
      HDivSurface, dense extract via `_to_dense(mat) = mat.COO()`.
    * `compute_centroids_areas_J(mesh, gf_J)` — per-tri J sampling
      via NGSolve `Integrate(..., element_wise=True)`, ready for the
      workpiece weak-coupling bridge.
* `src/radia/panels/calc_inductance.py:_solve_coil_bem_a` rewritten
  to call `compute_inductance_source_sink` + `compute_centroids_areas_J`.
  At maxh=0.005 the coil layer time dropped from 99.7 s to 2.3 s
  (43x faster) on the full weak-coupled IH pipeline.

### Files retired

* `src/radia/bem/efie_rwg.py` (-996 LOC; intree assembler + verification helpers)
* `tests/bem/test_coil_bem_a_dc_resistance.py` (-3 tests)
* `tests/bem/test_coil_bem_a_efie_golden.py` (-2 tests)
* `tests/bem/test_coil_bem_a_intree_match.py` (-3 tests)
* `tests/bem/test_coil_bem_a_order_convergence.py` (-2 tests)
* `tests/bem/test_coil_bem_a_rect_torus_loft.py` (-3 tests)
* `examples/coil_bem_a/` (entire directory — build123d demo, GMSH viz,
  PEEC-vs-BEM comparison scripts; the production path no longer needs
  these as showcase / orientation material)

### What was kept

* **Scalar BEM intree (Phase 1.9 C++ + 1.11 HACApK)** — kept and
  unchanged.  At N=595 intree scalar is 7x faster than ngsolve.bem
  with the COO fix, and supports HACApK ACA compression for large N.
  The HDiv asymmetry doesn't apply: H1 P1 has higher per-pair
  Galerkin product cost, which Phase 1.9's TaskManager parallel
  + admissibility cutoff offsets.
* **PEEC** path on `--coil-solver peec` — unchanged.

### Lesson learned

Benchmark new implementations against the alternative BEFORE shipping.
4.25.0 shipped intree BEM-A based on a single ngsolve.bem timing
(104 s @ N_J=5064) measured at a problem size 5x larger than the
panel uses.  At production N (300-1500), ngsolve.bem is ~0.4-6 s, not
slow.  See `memory/project_bem_a_ngsolve_chosen_2026_05_03.md`.

### Future Phase 6 (radia-mcp 0.39.x or radia 4.27.x)

Add curved-element support (`mesh.Curve(p)`) to the intree C++ scalar
BEM (`rad_bem_galerkin.cpp`).  ngsolve.bem already provides curved
geometry for the coil; only the workpiece (intree) needs the
extension.  Estimated 1.5-2 weeks: ~400 LOC C++ for isoparametric
Jacobian + Sauter-Schwab on curved triangles, ~50 LOC Python for
high-order node extraction, ~200 LOC tests.

## 4.25.1 — radia-heat console entry-point + Cubit-bypass launch documented

Released 2026-05-03.  Patch release.

### `radia-heat` console entry-point

The 4 standalone PySide6 panels (radia_ih, radia_em, radia_pcb,
radia_heat) are designed to run **without Cubit** — bring your own
`.vol` mesh (Cubit, Netgen-OCC, anywhere) and the panel runs
end-to-end.  Three of the four had `[project.scripts]` console
entry-points (`radia-ih`, `radia-em`, `radia-pcb`); `radia-heat` was
missing.  Now added: `pip install radia[gui]` registers all four
`radia-*.exe` launchers in `Scripts/`.

README "Standalone GUI Panels (no Cubit required)" section updated to
list all four launchers explicitly.

## 4.25.0 — BEM-A coil + unified inductance CLI + 6-method IH panel

Released 2026-05-02.  Surface-current coil solver lands as a peer of
PEEC; three scoped CLIs collapse into one.

### BEM-A coil solver (Phase C.1-C.5 shipped)

`radia.bem.efie_rwg` (~1 kLOC) — Weggler stabilized EFIE saddle on
HDivSurface RWG (= Lucy decomposition implicit), in-tree Sauter-Schwab
Duffy 4-D Galerkin assembler, source/sink port-driven self-inductance,
DC + AC SIBC closure for R extraction.  Tested against ngsolve.bem
oracle to 0.025-0.004 % on standard fixtures (5 test files, 13 cases).

Cross-method validation against PEEC golden on
`rect_torus_lofted_united.step` (8x6 mm rect, 50 kHz Cu):
- BEM-A converged 153.5 nH vs PEEC asymptote ~149 nH (n_peri → ∞)
- ~3.5 % residual gap is real modeling difference (BEM-A surface RWG
  resolves rect-corner current crowding; PEEC perimeter filaments do
  not).  Both methods are SIBC; both converge correctly within their
  discretisation classes.

### Unified `calc_inductance.py` (Refactor A1+B+C)

Single CLI replaces three scoped predecessors:

| OLD | NEW |
|-----|-----|
| `calc_peec_inductance.py` | `calc_inductance.py --coil-solver peec` |
| `calc_peec_bem.py` | `calc_inductance.py --coil-solver peec --vol …` |
| `calc_coil_bem_a_workpiece.py` | `calc_inductance.py --coil-solver bem-a --vol …` |

Dispatch on `--coil-solver {peec, bem-a}` × `--vol {present, absent}`
gives 4 modes (vacuum or weak-coupled, PEEC or BEM-A coil).  All
share the workpiece scalar BEM-SIBC block.

**Coupling terminology corrected**: previous "1-way forward" was
misleading because Telegen φ·(n·B) ΔL IS computed, capturing port-level
back-reaction even though coil J is fixed.  Renamed to **weak coupling**
throughout (CLI flag, panel labels, docstrings).  Strict one-way
(no ΔL) and strong coupling (coil J recomputed iteratively, FEM A-V)
remain distinct.

### IH panel: 6-method dropdown

`radia_ih.py` adds two BEM-A counterparts to the existing PEEC modes:

| Method | Coil | Workpiece | Notes |
|--------|------|-----------|-------|
| PEEC inductance (vacuum) | PEEC filament | — | existing |
| **BEM-A inductance (vacuum)** | BEM-A surface RWG | — | new |
| PEEC + BEM weak coupling | PEEC | scalar BEM-SIBC | renamed (was "1-way") |
| **BEM-A + BEM weak coupling** | BEM-A | scalar BEM-SIBC | new |
| PEEC coil + FEM Kelvin | PEEC | volumetric FEM-SIBC | unchanged |
| FEM A-V full | volumetric coil | volumetric FEM-SIBC | unchanged |

`coil_maxh` widget added (visible only for BEM-A modes).  19/19
`tests/panels/test_ih_panel_qt.py` PASS; 36/36 combined panel +
inductance + coil_topology golden tests PASS.

### PEEC unified OPEN/CLOSED cap-aware spine

New `radia.coil_topology` module (~280 LOC) — single source of truth
for OPEN vs CLOSED coil classification (cap detection, spine arc
parameters).  Consumed by:
- Path 2c (`_filaments_from_section_planes`) — `rect_torus_lofted_united`
  PEEC L 138.16 → 145.30 nH (filament was clipping ~14° off the 355°
  arc due to a `linspace(0, 2π)` fallback).
- Path 2b (`_filaments_from_circle_edges_per_station`) — single-turn
  coils now use topology spine; multi-turn helix guard
  (z_extent > 2 · median_r) keeps 3turncoil 422 nH on legacy NN-chain.
- Path 3 fallback — new `_centerline_from_topology_spine` helper
  inserted between `_centerline_from_torus_sweep` and
  `_centerline_from_open_spine`: closed-torus filament sweep
  178° → 358° (was tracing one half-arc seam edge).

Golden update: `peec_inductance_rect_united_50kHz_Cu.json`
138.16 → 145.30 nH, `regression_guard.min_L_nH` 120 → 130.

### Files

| Action | Path |
|--------|------|
| New | `src/radia/bem/efie_rwg.py` |
| New | `src/radia/coil_topology.py` |
| New | `src/radia/panels/calc_inductance.py` |
| Removed | `src/radia/panels/calc_peec_inductance.py` |
| Removed | `src/radia/panels/calc_peec_bem.py` |
| Removed | `src/radia/panels/calc_coil_bem_a_workpiece.py` |
| New | `tests/panels/test_inductance_golden.py` (6 tests) |
| New | `tests/test_coil_topology.py` (11 tests) |
| New | `tests/bem/test_coil_bem_a_*.py` (5 files, 13 tests) |
| Removed | `tests/panels/test_peec_inductance_golden.py` |
| Removed | `tests/panels/test_peec_bem_golden.py` |
| New | `examples/coil_bem_a/` (build123d demo + GMSH viz) |
| Updated | `src/radia/radia_ih.py` (6-method dropdown + dispatch) |
| Updated | `src/radia/panels/sync_registry.py` + `panel_registry.json` |

## 4.10.0 — Panel UX overhaul: font baseline, base-class unification, CoilBuilder wizards

Released 2026-04-26.  All-Python release; no C++ rebuild needed.
Focus is the panel runtime that EM, IH, and PCB share.

### Font baseline (readable on 2K)

The Qt OS default of 9pt Segoe UI was unreadable on the lab's 2K
displays (per `feedback_panel_vertical_space.md`).  After three turns
of bracketing with the user (9pt unreadable -> 11pt stingy -> 13pt
overshoot -> **12pt chosen**), `apply_panel_base_font(app)` is now
called from `run_app()` and from the test fixture; it inherits the OS
default family (Segoe UI on Windows, system sans-serif on Linux/macOS)
and bumps only the point size.

Hardcoded smaller fonts removed:
  - output text area: 9pt -> 12pt (Consolas)
  - status label (IH): "font-size: 11px" deleted -- de-emphasis is
    now colour-only (#888), font size inherits the baseline.

Regression guard: new `panel_qa.check_font_size_min` walks every
visible widget and fails if any renders below 10pt.  Width / height
thresholds rebudgeted (1400 RED / 1350 RED) for the larger glyphs.

### Panel base-class unification

EM and IH carried two identical 9-line copies of the section-header
helper.  One source of truth in `ModePanel`:

  - `_add_section(title, key=None)`           hoisted to base
  - `add_status_label(key=None)`              hoisted to base
  - `_method_combo` attribute convention      adopted by EM
                                              (was `_form_combo`)
  - `add_browse_action(key, label, callback)` new helper -- attach
                                              extra buttons to an
                                              existing browse row

PCB combo extraction `"0 (LU)".split()[0]` replaced with an explicit
`PCB_SOLVERS = {"LU": 0, "BiCGSTAB": 1, "HACApK": 2}` map.

### PEEC-inductance Window merged into IH

`radia_peec_inductance.py` (95-line wrapper) deleted.  The same
analysis is reached by selecting Method = "PEEC inductance (coil
only, STEP)" in the IH window.  IH's `__init__` auto-fills the STEP
field from the newest `*.step` / `*.stp` / `*.jou` in cwd when the
field is empty -- behaviour previously specific to the wrapper
window now applies to any IH launch in PEEC-inductance mode.

The Cubit launcher (`RadiaComp.cpp`) discovers panels by globbing
`radia_*.py` so the merged entry disappears from the menu without
any .ccl change.

### "New..." coil wizards (EM + IH)

Both panels expose a [New...] button next to the coil-source field:

  - **EM** -> writes a self-contained `.py` template (CoilBuilder
    racetrack with EDIT BLOCK at the top).  The user picks the save
    path, edits the 6 numbers (NI / START / WIDTH / HEIGHT /
    STRAIGHT / ARC_R), and points the panel's Coil script at it.
  - **IH** -> writes the same `.py` AND immediately runs it in-process
    to materialise a sibling `.step` (the format the IH PEEC modes
    consume).  User can later edit the `.py` and re-run
    `python coil.py` to refresh the `.step`.

Single COIL_TEMPLATE constant in `radia_gui_base.py` powers both.
Verified end-to-end: 84-segment closed racetrack, 2000 A, gap=1e-17 m,
.step = 125 KB OCC swept solid.

## 4.9.0 — EM canonical trio + 1/2 Kelvin Benchmark + MSC silent breakage fix

Released 2026-04-26.

### MSC silent breakage fix (production bug)

`calc_accel_msc.py` registered `add_material_args(include_custom=False)`,
so the EM panel's MSC mode (which sends `--material custom --mu-r <user>`
for the "mu_r (Linear)" material option) was rejected at argparse with
`error: argument --material: invalid choice: 'custom'`.  Clicking Run in
MSC mode produced no useful error in the GUI and no field result.  Fixed
by flipping to `include_custom=True` (matches the FEM Omega/A-Phi
convention in `calc_accel_magnet.py`).  Regression guard:
`tests/panels/test_em_msc_smoke.py` (6 sub-second static checks of the
MSC panel command vs argparse).

### Kelvin Benchmark — 1/2 sample joins 1/4

The EM panel's Kelvin Benchmark mode now ships **two** verification
samples (was 1/4 only):

| sample (frac) | reduction | error @ p=2 |
|---|---|---|
| `kelvin_benchmark_sphere_1_2.vol` | `{"y":"bn=0"}` | +1.07% |
| `kelvin_benchmark_sphere_1_4.vol` | `{"x":"bn=0","y":"bn=0"}` | +0.71% |

Build scripts refactored: parameterized core
`kelvin_benchmark_sphere_build.py` with `--frac 1_2|1_4|1_8` plus three
thin wrappers.  1/8 sphere build script ships for research only — the
.vol is **not** shipped because the magnetic-sphere-in-uniform-Hz BVP
fundamentally lacks 1/8 symmetry (z=0 mirror reverses Hzẑ source).
Realistic 1/8 EM workloads (C-yokes, dipoles) use the existing
`em/em_1-8_eighth.jou` C-yoke sample with the ELF `-x-y+z` convention.

Surprise (rho_min sweep): for compact-geometry Kelvin samples
(offset = 3R), the Periodic + sym BCs do most of the open-boundary
work — the `(R/r')^2` reluctivity is a small correction.  Capping
`Mu = mu_0` everywhere still gives 1/2 +0.34% / 1/4 -0.02% error.

### EM panel canonical trio (coil + .vol + BH)

Added two missing artifacts so the EM panel runs end-to-end on a
fresh checkout / pip install without requiring local Cubit
regeneration of every sample:

  - `src/radia/panels/samples/em_sample.vol` — half-z C-type dipole
    with auto-Kelvin (456k tets / 83k nodes).
  - `src/radia/panels/samples/em_sample_bh.txt` — 100-point steel
    BH curve (CEFC 2020 reference, 0..318 kA/m / 0..2.61 T).
    Full-precision source of the rounded built-in `STEEL_BH` table;
    functionally equivalent (relative diff < 2e-5).

Coil was already shipped as `em_sample_coil.py`.  Sub-second smoke
tests in `tests/panels/test_em_sample_artifacts.py` (8 tests) lock
ship-presence + parse-correctness for the trio.

### IH panel test refresh (4 stale test files repaired)

The 2026-04-19 IH panel restructure (4-method combo, separate
coil/wp material sections) left the test suite stale.  Rewrites:

  - `tests/panels/test_ih_panel_qt.py` — was 14/15 failing, now 17/17
    passing (method combo, per-method widget visibility, solver
    items, SIBC vs ESIM toggle, build_command roundtrip).
  - `tests/panels/test_panel_state_restore.py` — 6 tests rewritten
    from the obsolete 2-method labels to the canonical METHOD_*
    constants.  Round-trip widgets updated from retired
    `workpiece_mode` to surviving `impedance_model`.
  - `tests/panels/test_calc_inductance.py` — was failing at
    collection (ModuleNotFoundError on retired `bem_inductance`).
    `bem_inductance.py` now lives at
    `examples/induction_heating/bem_reference/`; tests `importorskip`
    cleanly when that dir is missing.
  - `tests/panels/test_ih_solvers.py` — same fix pattern as
    test_calc_inductance for the BEM precondition guard.

Also: `tests/panels/golden/em_eighth_mu1000.json` topology lock
refreshed for the post-2026-04-25 deterministic-anchor mesh
(ne 56289 → 56369, ndof 11695 → 11708; physics unchanged).
Slow physics regressions (`test_peec_bem_matches_2d_ref_coarse` 174s,
`test_fem_coilmesh_matches_2d_ref_gapped` 510s) both pass on the
shipped IH samples.

### Sample matrix documentation

Added `src/radia/panels/samples/README.md` indexing the canonical
samples per panel/method, identifying stale/research-only artifacts
that don't ship in CI-built wheels (gitignored .vol files).

### New panel logo

`src/radia/resources/2026_04_26.png` replaces the deleted
`radia_icon.png`.  `_icon_path()` now falls back to the most-recent
.ico/.png in `resources/` so future logo refreshes drop in without
a code change.

### radia-mcp 0.33.2

Updates `radia_ngsolve.kelvin_knowledge.benchmark_panel` topic to
document the 1/2 + 1/4 sphere samples + the rho_min sweep insight.
No API change.

## 4.7.0 — PEEC-inductance full geometry coverage + Japanese path + daemon speedup

Released 2026-04-22.

### PEEC-inductance — 4 coil topologies all STEP-only

The `radia_peec_inductance` panel mode ("PEEC inductance (coil only,
STEP)") can now extract a correct centerline from every coil topology
the lab uses, without requiring the user to supply an explicit `.jou`:

| Coil topology | Before 4.7.0 | After 4.7.0 |
|---|---|---|
| Circular gapped torus (single loop) | walker spine broken, L off 100× | ✅ TORUS analytical sweep, L = 85.10 nH (analytical 88.55, -3.9 %) |
| Rect 6 × 4 mm gapped torus | walker failed | ✅ CYLINDER analytical sweep, L = 88.15 nH |
| Multi-turn pancake loft (e.g. Kubota `3turncoil.stp`) | walker hung > 5 min, native crash | ✅ cross-section centroid NN chain with tangent continuity, L = 430.86 nH (matches `.jou` 426.25 nH to +1.1 %) in 12.9 s |
| `.jou` explicit centerline | | ✅ unchanged, fastest path |

New support code in `src/radia/coil_from_cad.py`:

- `_centerline_from_revolution_sweep()` — OCP `BRepAdaptor_Surface`
  extraction of TORUS / CYLINDER / CONE / REVOLUTION parameters,
  union of U-intervals = sweep angle, centroid of smallest PLANE
  face = spine radius.
- `_centerline_from_cross_sections()` — `GeomType.PLANE` face
  enumeration + dedupe near-duplicates + NN chain with tangent-
  continuity bias.
- `_auto_detect_cad_units()` — bbox-diagonal heuristic; handles
  Cubit STEPs that are sometimes mm-unit (values 10 – 1000) and
  sometimes m-unit (values 0.01 – 1) despite declaring
  `CONVERSION_BASED_UNIT('MILLIMETRE')`.

### Sibling `.jou` auto-preference (with format guard)

When `calc_peec_inductance.py --peec-step foo.step` is called and a
`foo.jou` sibling exists in the same directory, check whether the
`.jou` contains explicit `move Surface N x Y y Y z Z` lines before
switching to it.  A generic Cubit .jou (rect-sweep, brick setup)
no longer crashes the parser.

### Golden regression tests for PEEC-inductance

`tests/panels/test_peec_inductance_golden.py` — 3 tests locking in
the 3 canonical paths (TORUS analytical / cross-section centroid /
sibling .jou), each with a hard band + golden tolerance.

### IH panel workflows

PEEC+BEM (`calc_peec_bem.py`) and FEM A-V (`calc_fem_coilmesh.py`)
unchanged; verified on LAB + 100号機 that no regression.

### Packaging / deployment / documentation

- Cubit MCP daemon speedup (via `radia-mcp` 0.32.0 — see its CHANGELOG):
  cold start 30 – 60 s → 3 s.  VSCode restart 6 s → 0.01 s attach.
- `CLAUDE.md`: new Panel Samples Quality Policy + MCP Knowledge
  Placement Policy sections.
- `tests/panels/panel_qa.py`: 7 checks × every panel mode (deploy
  skill enforces).
- Deploy skill: per-mode Panel Mode Matrix + golden-range numeric
  assertions replace `status: "ok"` as the definition of "works".
- docs/research/IGA_CLN_DUAL_REDUCTION.md (LAB-local research plan).

### Known issues / non-goals

- Daemon auto-start across user logoff: Phase 2 deferred — Cubit
  GUI is desktop-bound so logoff-survival has no practical value.
- MCP server Python module hot-reload: still requires VSCode
  restart (OS-level limit, not addressable by daemon decoupling).

## 4.6.x — Panel QA + PEEC-inductance initial drop

(See git log before v4.7.0.  Adds the PEEC-inductance panel as a
new mode of the Radia-NGSolve launcher, plus Panel QA automated
layout checks.)
