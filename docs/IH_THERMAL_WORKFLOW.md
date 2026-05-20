# IH Thermal Workflow (radia 4.59.0+)

Complete EM → Thermal pipeline in a single `radia-ih` panel, from a
PEEC / FEM-Kelvin coil + workpiece solve all the way to the
workpiece temperature distribution.

## Phase A: EM solve → q_surf .sol

Pick a method in the **radia-ih** panel that emits a surface heat
flux .sol on the workpiece SIBC surface:

| radia-ih method | Emits q_surf? | Driver script |
|---|---|---|
| PEEC inductance (coil only) | NO | calc_inductance.py |
| BEM-A inductance (coil only) | NO | calc_inductance.py |
| PEEC + BEM weak coupling | **YES** | calc_inductance.py --coil-solver peec |
| BEM-A + BEM weak coupling | **YES** | calc_inductance.py --coil-solver bem-a |
| PEEC coil + FEM wp (SIBC) + Kelvin | **YES** | calc_fem_kelvin.py |
| Full simulation (FEM A-V + wp SIBC + Kelvin) | **YES** | calc_fem_coilmesh.py |

The `qsurf_sol` JSON key in the result holds the absolute path to
`<stem>_qsurf.sol`.  The matching mesh is at `<stem>_fem.vol`.

## Phase B: Thermal solve

Pick **Method → "Thermal (heat transfer from saved q_surf .sol)"**
in the same `radia-ih` panel.  The EM-side sections (Drive, Coil
material, Coil geometry, Workpiece material, Workpiece impedance,
Linear solver, Advanced) hide; the embedded HeatPanel becomes the
active control surface.

### Two ways to enter the Thermal method

1. **Manual**: pick "Thermal..." from the Method dropdown directly.
   Fill the qsurf .sol + EM .vol fields by hand.
2. **Chain shortcut**: after a successful EM solve produced
   `qsurf.sol`, click **"Run thermal..."** in the action row.  The
   button switches the dropdown to Thermal and pre-fills
   `qsurf_sol` + `em_vol` from the JSON result.  Set heat source =
   "Spatial q_surf .sol (from IH)" automatically.

## The `.sol` + `.vol` contract (v4.58.0 strict)

Both files MUST be supplied; the panel's Run button stays disabled
until both exist on disk.

* **`qsurf_sol`** -- NGSolve `.sol` written by `gf_q.Save(...)` in
  the EM solver.  Pure binary coefficient vector; no embedded mesh,
  no FES-order header.  Useless without the matching mesh.
* **`em_vol`** -- the EM `.vol` the `gf_q` GridFunction was saved
  against.  The thermal solver rebuilds
  `H1(em_mesh, order=qsurf_order)` and `Load`s the .sol into it.
* **`qsurf_order` MUST equal the EM solve's `fes_order`**.  No
  metadata exists to auto-detect this; a mismatch reads garbage
  silently.  Default both = 1 (the radia panel default).  When the
  EM solve uses `--fes-order 2`, set `--qsurf-order 2` in the
  Thermal sub-panel.

How the contract is enforced:

| Layer | Behavior |
|---|---|
| `calc_heat.py --qsurf-sol PATH` without `--em-vol` | Raises `ValueError` with a hint at the typical sibling .vol path |
| `radia._heat_panel.HeatPanel.is_runnable()` | Returns False unless both .sol AND .vol exist on disk |
| `radia_ih._on_run_thermal` chain helper | Tries 3 strategies (msh_file stem, qsurf stem `_qsurf.sol`→`_fem.vol` swap, JSON `wp_vol`) to auto-locate the EM .vol; emits a hint line when all 3 fail |

Auto-locate from the .sol stem was **removed in v4.58.0** per the
"No Fallbacks" policy — silent fallbacks here had been masking
mesh-vs-coefficient mismatches.

## Mesh: 3D volume vs 2D axisymmetric

The HeatPanel's "Mesh type" combo routes to one of two solvers:

* **3D volume** → `calc_heat.py`.  Arbitrary 3D workpiece.
* **2D axisymmetric (r, z)** → `calc_heat_axisym.py`.  10-100×
  faster on rotationally symmetric workpieces (cylinder, stepped
  shaft).  Cross-mesh q_surf transfer is φ-averaged so a slightly
  non-axisymmetric coil (gapped torus) still produces a physically
  sensible q.

Both solvers consume the same `qsurf_sol` + `em_vol` pair.

## Heat source mode

* **Uniform q_surf [W/m²]** -- a single scalar applied across the
  whole heating face.  Cheap, useful for first-cut feasibility runs
  and smoke testing.  Rotation has no effect (constant in space);
  the panel warns when `rotation_rpm > 0` is combined with this.
* **Spatial q_surf .sol (from IH)** -- load the .sol that the EM
  solve emitted.  Preserves hotspot distribution.  Rotation
  re-samples q_surf on the body frame each timestep.

## Workpiece rotation (v4.58.0+)

Set `Rotation [rpm]` > 0 in the Thermal sub-panel.

* **3D solver** (`calc_heat.py`): the workpiece body spins around
  the z axis.  q_surf is re-projected each timestep at the body's
  instantaneous angle.  Mesh / FES / mass / stiffness remain fixed;
  only the LinearForm RHS reassembles (which it already did for
  the convection term, so per-step overhead = one re-projection ≈
  10 ms on typical meshes).  Default `max_iter=60`, `tol=1e-3`,
  `relax=0.3`.
* **2D axisymmetric solver** (`calc_heat_axisym.py`): rotation is
  implicit in the axisymmetric assumption.  The rpm value is
  recorded as metadata; no time-loop change.

The body-frame projection at angle θ:

```python
c, s = math.cos(theta), math.sin(theta)
for vnr, (xb, yb, zb) in zip(surf_vnrs, surf_xyz):
    xw = xb*c - yb*s
    yw = xb*s + yb*c
    em_mip = em_mesh(xw, yw, zb)        # world-frame lookup
    val = gf_q_em(em_mip)               # q_em evaluated there
    gf_wp_q.vec.FV()[vnr] = float(val)
```

Validated against synthetic 2-wire setups (q_em(x,y,z)=x at
theta=0/π/π/2 reads body-(1,0,0) as +1/-1/0).  See
`tests/panels/test_heat_rotation.py`.

## Output

The 3D solver writes:

* JSON to stdout with `T_max_C`, `T_min_C`, probe history, time
  history, `Q_input_J`, runtime stats, AND the saved file paths
  (`T_sol_file`, `heat_vol_file`, `msh_file`, `vtu_files`).
* **T `.sol` + companion `.vol`** -- the final temperature
  GridFunction is ALWAYS saved as a NGSolve `.sol` next to the
  workpiece thermal mesh.  Symmetric with the EM side's qsurf.sol
  contract:
    - With `--msh-output`: writes `<msh-stem>_T.sol` +
      `<msh-stem>_heat.vol` (a fresh companion mesh) alongside
      the GMSH `.msh`.
    - Without `--msh-output`: writes `<wp-stem>_heat_T.sol` next
      to `wp.vol` (re-use wp.vol itself as the companion mesh).
  Both paths are reported in the JSON as `T_sol_file` and
  `heat_vol_file`; the IH Summary lines them out next to GMSH /
  VTU paths so the user sees them at a glance.
* GMSH `.msh v4.1` (per the lab standard) via
  `gmsh_post_export.vol2msh`.  Bundled fields: `T_C` (volume
  scalar, per-vertex), `q_surf` (surface scalar on the heating
  face).
* Per-timestep `.vtu` files when `--vtu-prefix` is supplied
  (transient animation).
* Optional `.csv` probe history when `--probe-point x,y,z` +
  `--csv-output` are both set.

### Reloading the T `.sol` for later evaluation

```python
from ngsolve import Mesh, H1, GridFunction
wp_mesh = Mesh("workpiece_thermal.vol")    # or <stem>_heat.vol when --msh-output was used
fes_T   = H1(wp_mesh, order=1)             # MUST match the solve's --fes-order
gfT     = GridFunction(fes_T)
gfT.Load("workpiece_thermal_heat_T.sol")   # or <msh-stem>_T.sol
# Now gfT is the final temperature field on the thermal mesh:
T_at_point = float(gfT(wp_mesh(0.0, 0.0, 0.005)))   # sample at body point
```

Same contract as qsurf.sol: the .sol is a coefficient vector only;
the matching `.vol` + matching FES order are required to reconstruct
the GridFunction.

## Trouble shooting

| Symptom | Cause | Fix |
|---|---|---|
| Run button disabled in Thermal section | Either `wp_vol` is empty, or `qsurf_sol` empty in spatial mode, or `em_vol` empty in spatial mode | Browse to all three files. |
| `--em-vol is required` error | Manual CLI invocation of `calc_heat.py --qsurf-sol PATH` without the companion .vol | Pass `--em-vol <stem>_fem.vol` explicitly |
| `T_max ≈ T_initial` (no heating) | wp surface vertices fell outside the EM mesh → `n_fail` count in the log | Verify the EM and thermal meshes describe the same physical workpiece geometry (check coordinate origin + units) |
| `T_max` jumps when switching `fes_order` to 2 | `qsurf_order` was not updated to match | Set `qsurf_order = 2` in the Thermal sub-panel to match the EM solve's basis order |
| Rotation has no effect on `T(t)` | Heat source = Uniform (constant in space) OR `rotation_rpm` left at 0 | Set heat source = Spatial AND `rotation_rpm > 0` |

## Removed: `radia_heat` standalone (4.62.0+)

The standalone `radia_heat.py` module and the `radia-heat` console
script were removed in radia 4.62.0.  Heat analysis is the
"Thermal" Method choice in `radia-ih`; the HeatPanel sub-widget
now lives at `radia._heat_panel` as an internal implementation
detail of the IH panel.  The pre-4.59 standalone window is gone
and there is no public CLI replacement for `radia-heat` -- launch
`radia-ih` instead, then pick Method = "Thermal".

Pre-4.59.0 shortcut behavior: in 4.59.0-4.61.0 ``radia_heat.py
main()`` was a deprecation stub that redirected to ``radia-ih``.
That stub is gone in 4.62.0; calls to `radia-heat` from old
scripts will fail with "No module named 'radia.radia_heat'"
(import) or "command not found" (CLI).  Update shortcuts to
launch ``radia-ih``.

Programmatic imports: replace
``from radia.radia_heat import HeatPanel, HEAT_SRC_SPATIAL`` with
``from radia._heat_panel import HeatPanel, HEAT_SRC_SPATIAL``.

## Test coverage

* `tests/panels/test_heat_rotation.py` -- unit tests for the
  rotation projection math on a synthetic unit cube
  (`q_em(x,y,z)=x` reads +1/-1/0 at body-(1,0,0) for θ=0/π/π/2).
* `tests/panels/test_heat_chain_golden.py` -- e2e chain test
  invoking `calc_heat` via subprocess.  Slow-marked.
* `tests/panels/panel_qa.py` registry includes `ih_thermal` which
  exercises the integrated Thermal method through `IHWindow`
  (rendering / font / layout checks).

## See also

* [`docs/peec/VOLUME_PEEC_DESIGN.md`](peec/VOLUME_PEEC_DESIGN.md)
  -- the deferred Volume PEEC design (radial filaments to capture
  proximity effects beyond perimeter PEEC).
* [`docs/esim/R_MISMATCH_PEEC_VS_BEMA.md`](esim/R_MISMATCH_PEEC_VS_BEMA.md)
  -- the perimeter-PEEC vs BEM-A R discrepancy with proximity
  factor 1.22 ceiling at 150 kHz.
* radia-mcp tools: `ih.rotating` (calc_heat.py rotation impl),
  `radia_ngsolve.ngsolve` Section 18c (qsurf-style H1 GF
  save/load).
