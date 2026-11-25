"""PEEC route knowledge for induction heating (radia_mcp.ih).

PEEC and ESIM are the main Radia routes of the IH pack. Independent FEM
cross-checks of these routes use the FEM knowledge of ngsolve-mcp; general
PEEC formulation lives in the ``peec`` pack.
"""

from radia_mcp.matrix_solvers.direct_solvers_knowledge import SOLVER_CAPACITY


IH_PEEC_ROUTES = """
# IH PEEC routes

PEEC (coil filaments) and ESIM (nonlinear surface impedance) are the main
Radia routes for induction heating. The PEEC coil enters the IH workflow in
two production couplings:

- PEEC + BEM-SIBC (weak coupling): coil-only or uniform-SIBC workpiece,
  `python -m radia.panels.calc_inductance --coil-solver peec --vol <wp.vol>`.
  Yields L (Telegen delta-L) and workpiece loss P_wp. See topic `bem_sibc`.
- PEEC + FEM workpiece + Kelvin: prescribed PEEC filament source driving an
  HCurl FEM-SIBC workpiece, `python -m radia.panels.calc_fem_kelvin`.

Nonlinear steel or ferrite workpieces use ESIM with Karl iteration
(`ih_esim`). Linear SIBC architecture and solver boundaries are in `ih_sibc`.
General PEEC formulation, bundles and AC loss are in the `peec` pack.

Independent check: compare against a FEM route (volume coil A-V with
resolved skin, `induction_heating(topic="av_coil_sigma")`, and the FEM
knowledge of ngsolve-mcp). A PEEC result is not validated by its own
approximation.
"""

IH_PEEC_BEM_SIBC = """
# PEEC (coil filament) + BEM-SIBC (workpiece) — production path for rotating WP

## POLICY: BEM-A = impedance-EFIE, solved by COCR or HACApK-COCR

**POLICY (2026-07-03, Sugahara):** the BEM-A coil solver
(``--coil-solver bem-a``, ``radia.bem.coil_inductance_ngsolve``) uses the
**impedance-EFIE** formulation and ONLY that, and it is solved by
**COCR or HACApK-COCR** -- the two ``--coil-saddle-solver`` choices:

- ``cocr``        : div-free loop reduction + complex-symmetric COCR, dense
                    SL matvec (``auto`` picks this for large systems).
- ``hacapk_cocr`` : the same COCR with the HACApKBEMManager-compressed
                    O(N log N) SL matvec (accuracy ~3e-7; identical R/L).
                    Since 2026-09-21 it accepts any ``--coil-fes-order``: the
                    cluster tree takes one point per DOF (edge midpoints for
                    edge DOFs, centroids for face DOFs) instead of assuming
                    RT0.  Raise the order to p-converge current crowding on
                    strongly curved conductors; on the beak-fin fixture one
                    order step moved R by +1.69% where curving the geometry
                    moved it by +0.067%, so ``--maxh`` alone does not reach
                    a rounded tip whose element size is set by its curvature.

(Dense ``lu`` stays the small-N direct option; ``gmres`` / ``minres`` are kept
only for comparison -- unpreconditioned GMRES stalls on the indefinite AC
saddle.)  The Leontovich surface impedance ``Z_s = (1+j)/(sigma*delta)`` sits
INSIDE the saddle system::

    [ jw*mu0*SL + Z_s*M   D^T ] [J]   [0]
    [ D                    0  ] [p] = [g]      (complex, omega > 0)

so the recovered J is the FINITE-IMPEDANCE surface current, and
``R = Re(Z_s) * (J^H M J)`` (Leontovich SIBC dissipation),
``L = mu_0 * (J^H SL J)`` (EXTERNAL inductance -- the internal surface
reactance ``Im(Z_s)*(J^H M J)/omega`` is deliberately NOT folded into L).
At ``omega == 0`` it reduces to the real vacuum-L saddle (R = 0).

**The historical PEC post-hoc formulation is REMOVED and must not return.**
It solved a perfect-conductor saddle then integrated ``R = Re(Z_s)*|J|^2 dS``
afterwards; the PEC J concentrates singularly at near-contact turn gaps /
edges where it varies below the skin depth and the Leontovich integral
breaks down, over-estimating R ~3x on tightly-wound coils (example 3-turn
pancake: PEC 15.14 mOhm vs the physical 4.63, the latter confirmed by
volume PEEC 3.7 / perimeter PEEC 4.5 / analytic proximity 4.8).  On smooth
geometry both agreed (isolated straight wire = closed-form Bessel to <1%),
so the removal loses nothing.

**Scalable solver = COCR** (two ``--coil-saddle-solver`` choices: ``cocr``
[dense SL matvec] and ``hacapk_cocr`` [HACApK-compressed SL matvec]; ``auto``
selects ``cocr`` for large systems: > 5000 tris AC / 10000 DC).  It reduces
the saddle to the divergence-free (loop / stream-function) subspace via the
sparse projector onto ker(D), giving the complex-symmetric operator
``Pi A11 Pi`` that COCR (Sogabe-Zhang, UNCONJUGATED inner products) solves in
~24 MESH-INDEPENDENT iterations -- EXACT vs the dense LU (dR/dL < 0.01% on the
gapped torus) and replacing BOTH the O(N^3) LU and the unpreconditioned GMRES
that stalled (~1e5 matvecs) on the indefinite AC saddle.  ``hacapk_cocr`` uses
the HACApKBEMManager-compressed O(N log N) SL matvec (accuracy ~3e-7, the
bem_sibc_solver pattern); ``cocr`` and ``hacapk_cocr`` give identical R/L.
COCR does NOT fit the raw saddle
(structural breakdown -- rhs = [0; g] lives in the constraint block so the
initial r^T A r = 0) nor the Schur complement (diverges); the div-free
reduction is what makes it fit.

## L differs from PEEC by ~4% -- BENIGN, and curve order does NOT close it

R is the accepted metric and it AGREES (impedance-EFIE 4.63 mOhm ~= perimeter
PEEC 4.5 / analytic 4.8, all in the ~4.5-5 band).  The residual L gap --
BEM-A EXTERNAL L 411.6 nH vs PEEC 430 nH (~4.3%) on the example 3-turn coil --
is a KNOWN, BENIGN geometry/convention difference, LEFT AS-IS (Sugahara
2026-07-03: signed off because R agrees).  It is NOT a defect, and raising the
mesh curve order does NOT close it:

- MEASURED (torus 150 kHz, fixed RT0 DOF n_J=5745): curve_order 1/2/3 give
  L=87.2702 nH and R=1.04531 mOhm BIT-IDENTICAL.  ``mesh.Curve(p)`` is a NO-OP
  in this path because the in-memory surface extraction
  (``_extract_surface_mesh_filtered``) strips the CAD/curvedelements
  association, so there is no geometry to project the mid-side nodes onto.
  Raising geometry order would require a Cubit ``export netgen ... order N``
  surface ``.vol`` loaded directly.
  UPDATE 2026-09-21: the coil path no longer extracts.  The volume ``.vol``
  is passed whole and ``compute_inductance_source_sink`` compresses the
  space, so the export-time curving order (and ``mesh.Curve(p)`` on an
  in-memory OCC mesh) now reaches the solve.  Measured on the beak-fin
  fixture: geometry order 1 -> 2 moved R by +0.067%, confirming the
  conclusion above; the lever that does move a curved conductor is the
  BASIS order (``--coil-fes-order``: one step +1.69%).
- Even properly applied, curve order is NOT the lever: EXTERNAL L is
  loop-dominated (turn radius, enclosed area, N^2 -- already resolved by flat
  elements), and round-wire faceting enters L only logarithmically
  (ln(8R/a)) = a sub-1% effect.  The ~4% gap is the "surface mesh != PEEC
  filament geometry" mismatch + PEEC finite n_peri + external-L(BEM) vs
  bundle-L(PEEC) convention -- none touched by curve order.
- 430 nH is NOT validated ground truth (it is PEEC's own approximation; the
  FEM A-V 3-D reference is still TBD).  If L ever needs reconciling: match the
  SAME STEP geometry, h-refine the surface mesh, then anchor to FEM A-V
  (``calc_fem_coilmesh.py``) or measurement -- NOT curve order.

Refs: ``docs/peec/VOLUME_PEEC_DESIGN.md``,
``docs/esim/R_MISMATCH_PEEC_VS_BEMA.md``,
``validation_test/bem/test_coil_bem_a_impedance_efie.py``.

## SIBC reaction, heating and genus-1 workpieces (2026-09-15)

The verified weak workpiece route uses the SOLVED total surface magnetic
field for heating. The former incident Biot-Savart pattern normalized to
global loss has been removed. NGSolve assembles a positive lumped P1
projection of q = 0.5 Re(Zs) |H_t|^2. Its integral must equal the BEM loss;
an incomplete parent-mesh mapping or failed qsurf export raises.

Coil reaction uses COMPLETE SIBC reciprocity: the magnetic phi.B term
PLUS the electric boundary term. For a single-valued surface potential,
the added inductance is Zs/(i omega I^2) psi.T K phi (no conjugation).
It is not obtained by defining resistance from dissipated power.
Outputs identify ``telegen_form=complete-SIBC-reciprocity``,
``wp_reaction_power_W``, ``wp_power_balance_relative_error``,
``wp_power_balance_tolerance``, and
``qsurf_method=solved-total-field-lumped-P1``.
The runtime refuses negative reaction power, nonfinite results or a
power imbalance above 10%; this safety ceiling is NOT an accuracy claim.
Numerical acceptance uses tighter mesh-converged FEM/BEM comparisons.

Supported weak postprocessing is surface P1. Genus-0 weak coupling accepts
uniform Zs or per-BND-triangle ESIM with coefficient-weighted stiffness.
This restriction is on the EM BEM basis, NOT the thermal FEM order;
standard NGSolve H1 order 2 remains the thermal route.
Higher-order BEM remains rejected. Local loss and complete reaction use the
same source-weighted form; heat is conservatively projected to P1. No fallback is selected.
Scalar ESIM remains subject to its constitutive convergence checks;
linear reference tests alone do not certify arbitrary nonlinear BH data.

For genus 1, ``--wp-loop-dof auto`` activates the SIBC or per-panel ESIM,
intree-dense, P1 extension. Both spatial heat and reaction include the
harmonic carrier, not just the reported total loss. The reciprocal
current.A contribution on the cut-open field retains the cut term.
The loop closure uses distributed Galerkin electric work and magnetic
work including both tangential current and normal flux. Magnetic work
acts on the scattered field; electric work and heat use the total field.
Work pairings are reciprocal; the mixed BIE matrix is not symmetric.
This path requires flat undeformed surface triangles. The additional P0
single layer is cached per geometry. ``--wp-loop-work-backend auto`` selects
native Galerkin FMM at >=512 faces and native direct products below that;
``dense`` and ``fmm`` explicitly select either route in weak and both strong
workflows. Separated interactions are compressed; direct near/singular
quadrature and the complete closure are unchanged. The result's
``wp_loop_work_diagnostics`` records the route, FMM controls and NGSolve
version. Unavailable FMM raises with the dense option, without fallback.
Native FMM storage is unknown, never zero; use measured process peaks.
The P1 BIE stays dense and the 14000-face loop guard remains. Compression
heat-change acceptance is <=0.1% at the same mesh, not a physical-accuracy
bound. Carrier, seam, quadrature and surface refinement still require checks.
Unsupported weak handle combinations (scalar ESIM, HACApK, P2, genus >= 2)
raise before BEM assembly. There is no ``off`` fallback. Surface winding
is established by face-BFS plus signed volume, never centroid flipping
of an inner wall. The P1 incident-potential surface-Poisson residual
must stay below 10%.

The repository-owned acceptance is
``validation_test/eddy_current_analytical_validation/validate_ih_sibc_closure.py``:
axisymmetric FEM-SIBC versus BEM on solid and bored cylinders, separate
source power and dissipation, local complex magnetic field, surface
refinement, and FEM exterior-size convergence. The full-conductor FEM
power baseline and analytic shorted-ring test are complementary checks.
No external case files or reported ratios are acceptance inputs.

Scalar strong coupled coil solvers use complete mutual reaction from
TOTAL body surface current and the discrete incident electric projection.
The genus-1 mode is solved on every coil-current iterate; post-convergence
loop correction is no longer used. The surface-coil single-layer matrix
is unscaled, so its reaction RHS includes 1/mu0. Coil fixed-point and body
residuals certify the same returned state. Body heat is checked against
complete body reaction; port power includes coil self loss, and port loss
increment includes the CHANGE in coil self loss. It is not inferred from
body heat alone. The surface-coil strong model is ideal (zero coil ohmic
loss); the filament route uses its R_f/Zs_fil loss. Each air baseline uses
the same model as its coupled result. Fixed panel Zs is supported by the Python/CLI strong surface-coil and
filament-coil routes for flat P1 dense genus-0/1 bodies. Values are used
face by face; the area-weighted mean is reporting metadata only. Nonlinear
per-panel ESIM connects both strong coil routes through a damped material
Picard loop at the actual terminal current. Inner coupling tolerance is
<= esim_tol/10. The final complete EM re-solve and direct ALL-panel cell
certificate share body/coil/loop and reaction/heat gates. Port quantities are
secant-at-accepted-current; changing current needs a new material solve.
Outer and inner iteration histories are separate. The native Simulink strong
thermal handoff remains guarded.
Prescribed-source genus-0 phi.B quadrature is unchanged: it converges to
Js.A work but is not identical on a finite mesh. Strong coupling has its
own self-authored circuit and surface checks; weak acceptance and native
Simulink certification do not transfer to it.
Per-panel ESIM defaults to `--esim-panel-evaluator table` for weak genus-0/1
and strong BEM/PEEC. The adaptive log-H table reduces outer cell evaluations;
all-panel direct final certification is unchanged. Select
`--esim-panel-evaluator direct` for the reference evaluator. Budget, refinement
or field-range table failures raise with this guidance; no silent fallback.
The selected mode is recorded in `esim_panel_evaluation.mode`.


## When to use this

- Workpiece moves/rotates relative to coil (different relative position per step)
- Coil geometry is fixed
- Linear or nonlinear SIBC on workpiece (delta << WP size)
- Frequency range where filament model captures coil AC redistribution

## Advantages

- **No coil volume mesh needed** — coil is a filament bundle from
  `CoilBuilder.to_filaments()` or `coil_from_build123d_sweep()`.
- Only workpiece surface mesh + per-step rebuild as WP moves.
- BEM solver construction cached once per WP mesh (`IHWorkpieceContext`).
- L from PEEC `Z_port` is mesh-independent (matches analytical to <1%).
- Per-coil-variant eval ~3 seconds (MNA + Biot-Savart + BEM linear solve).

## Pipeline (validated 2026-04-18)

```python
from coil_builder import CoilBuilder
from peec_bundle import build_bundle_solver
from radia.bem_sibc_solver import (ScalarBIESIBCSolver,
                              compute_phi_inc_from_filaments)

# 1. Coil geometry via CAD-first API (NO mesh)
coil = (CoilBuilder(current=I_total)
        .set_start([R_major, 0, 0])
        .set_cross_section(width=SIDE, height=SIDE)
        .add_arc(radius=R_major, arc_angle=355))
paths, _ = coil.to_filaments(nw=3, nh=3, frequency=f,
                              sigma=5.8e7, mu_r=1.0, n_arc=40)
info = coil._last_filament_info
dw, dh = info['w']/3, info['h']/3

# 2. PEEC MNA -> per-filament complex AC currents
solver_peec, seg_of_fil, _, _ = build_bundle_solver(paths, dw, dh, 5.8e7)
I_branch = solver_peec.compute_branch_currents(f, [I_total+0j])
I_peec = [np.mean([I_branch[s] for s in segs]) for segs in seg_of_fil]

# 3. Inductance from port impedance (MESH-INDEPENDENT)
Z_port = solver_peec.compute_port_impedance(f)
L = Z_port.imag / (2*pi*f)
R_ac = Z_port.real

# 4. phi_inc at WP surface via Biot-Savart line integrals
wp_pts = np.array([wp_mesh.vertices[v].point for v in range(wp_mesh.nv)])
phi_inc = compute_phi_inc_from_filaments(wp_pts, paths, I_peec, n_quad=20)

# 5. BEM-SIBC solve on workpiece (CACHE THE SOLVER)
bem = ScalarBIESIBCSolver(wp_mesh)  # one-time ~10s
sol = bem.solve(phi_inc, Z_s, omega)
H_t_rms = sol['H_t_rms']
P_total = sol['P_density'] * sol['area']
```

## Validation vs A-V FEM truth (2026-04-18)

Geometry: gapped torus R=30mm, a=3mm, 5 deg gap. Steel WP R=25mm H=25mm.
Frequency 7 kHz (Cu delta = 0.79 mm, delta/a = 0.26 — weak skin effect).

| Method | L [nH] | H_t [A/m] | P [W] |
|--------|--------|-----------|-------|
| A-V FEM (coil-sigma) — truth | 100.25 | 10.36 | 4.95e-4 |
| PEEC + BEM-SIBC | 88.57 | 10.67 | 5.24e-4 |
| Agreement | -11.6% | **+3.0%** | +4.4% |

L difference reflects FEM truncation vs analytical; H_t/P agree within method
differences (FEM volume vs BEM surface on different WP mesh).

## mu_r support

Linear SIBC in the workpiece: `Z_s = (1+j)/(sigma*delta)` with
`delta = sqrt(2/(omega*mu_0*mu_r*sigma))`. Pass mu_r via the material setup
in `IHWorkpieceContext(mesh, freq, sigma, mu_r=100, ...)`.

Nonlinear BH cell problems can produce spatially varying impedance, but the
current ScalarBIESIBCSolver accepts uniform Z_s or PanelSurfaceImpedance on
surface P1. Untagged nonuniform arrays fail loudly: vertex values are not panel
values. Source-side weighted stiffness K_gamma and local loss integration use
one Zs per BND triangle; --esim-per-panel connects this to the genus-0 weak
coil/heat route, weak genus-1 and strong genus-0/1 dense flat P1 EM routes.
Genus-1 and strong nonlinear iterations use damped Picard with simultaneous
final all-panel direct constitutive certification. This implementation does not certify arbitrary magnetic materials.

## Coil side: copper (linear, mu_r=1) is the default assumption

PEEC filament MNA uses sigma_coil for AC R and M. Copper at typical IH
frequencies is non-magnetic (mu_r=1). For ferromagnetic coil materials,
PEEC filament model would need revision (not currently supported).

## Panel integration: calc_inductance.py (v4.25.0+ unified Layer 4)

The `src/radia/panels/calc_inductance.py` script is the Layer-4
subprocess for the IH panel's "PEEC + BEM weak coupling" method.  It
unified three legacy scripts (`calc_peec_inductance.py`,
`calc_peec_bem.py`, `calc_coil_bem_a_workpiece.py`) into one CLI with a
``--coil-solver {peec,bem-a}`` switch.  Invocation:

    python calc_inductance.py \\
        --coil-solver peec --coil-step <coil.step> \\
        --vol <wp.vol> --wp-label sibc \\
        --sigma <S/m> --mu-r <...> --half-thickness <m> \\
        --frequency <Hz> --current <A>

It wraps the demo pipeline above into a .vol-based / JSON-stdout
interface.

### BND extraction from wp-as-hole geometry (CRITICAL)

When the workpiece is modeled as a HOLE in the air volume with a sibc
sideset (closed-torus / skin sample convention), multiple FaceDescriptors
may share the same 'sibc' name (e.g. air_top + air_bot halves after
webcut).  Each FD's BND element winding gives an outward normal relative
to its own domin, which CAN BE OPPOSITE between FDs.  A blanket flip is
WRONG — 77% of triangles needed flipping in the closed-torus sample,
23% did not.

Robust fix (implemented in `_extract_bnd_only_inline`, formerly in the
deleted `calc_peec_bem.py`):
1. Collect the wp-surface vertex cloud.
2. Compute wp_centroid = mean(coords).
3. For each triangle, check sign of dot(normal, centroid_to_face).
4. Flip winding if inward.

Without this, BEM P output is silently wrong by 2-20x (not a failed
solve — a plausible-looking but incorrect number).

### Golden regression test (policy)

`validation_test/panels/test_inductance_golden.py` runs `calc_inductance.py` on
`ih_peec_bem_coarse.vol` at 7 kHz Cu and asserts:
- L_coil within 1% of captured 88.57 nH (geometry-only, tight)
- P_wp within 15% of 2D SIBC reference 6.63e-5 W (mesh-sensitive)

Expected values in `validation_test/panels/golden/peec_bem_coarse_7kHz_Cu.json`.
Re-generate the sample via Cubit headless:
`coreform_cubit -batch -nojournal -nographics -input
src/radia/panels/samples/ih_peec_bem_coarse.jou`.

### Sibling: calc_fem_coilmesh.py (2026-04-19 FINAL, A-V + wp SIBC)

Gapped-torus A-V formulation following the MCP-documented
`INDUCTION_HEATING_AV_COIL_SIGMA` pattern.  Coil is volumetrically
meshed with sigma; phi (scalar potential) is defined on coil material
only, Dirichlet 1 on 'source' / 0 on 'sink' port sidesets.  Volume-
integral current extraction `I_out = int J . grad(psi_n) dV`.  WP
SIBC Robin as usual, Kelvin exterior.

**Design journey that got us here (retired options)**:
- Neumann-K on coil surface alone (2026-04-19 am): L collapsed 10x
  because Robin SIBC cancels K_imposed exactly.  **Retired.**
- Neumann-K without Robin (pure Dirichlet-source): got P_wp within
  2D SIBC ref, but P_wp matched PEEC+BEM only to 24% — the uniform-K
  surface source misses coil proximity effects.  **Retired.**
- Dowell thin-wire analytical P_coil: useful as a validation reference, but not
  a production-block primitive because it hides mesh dependence. The retained
  reference belongs under `validation_test/induction_heating/`.

### Validation (2026-04-19, ih_fem_kelvin_skin_fine, Cu 7 kHz)

| Method | P_wp [W] | L [nH] | P_coil [W] |
|--------|----------|--------|------------|
| 2D axisym (ref) | 6.63e-5 (SIBC) | 50.54 | ~8.4e-5 (implied) |
| PEEC+BEM | 6.48e-5 (-2.2%) | — | — |
| FEM A-V   | **6.54e-5** (-1.3%) | 47.32 (-6.4%) | 1.48e-4 (+76% mesh-sens.) |

**P_wp agreement PEEC+BEM vs FEM A-V: 1%** (!) — the primary IH
engineering metric.  L is within 10% (coil mesh-sensitive).  P_coil
is volumetric `|J|^2/sigma` and mesh-sensitive; captured value with
tolerance 15% in golden. Use the checked Dowell analytical reference under
`validation_test/induction_heating/` when a coil-loss comparison is required.

### .vol / .jou requirements

- 'coil' material, VOLUMETRIC-MESHED with `volume size <= delta/3`
  (ideally; skin_fine sample uses 0.16mm surface + 0.5mm volume for
  Cu 7kHz delta=0.79mm).
- 'source', 'sink' sidesets on the two gap-face circles.  Real IH
  coils have physical port terminations; use those.  Closed-torus
  topology cannot drive A-V.
- 'sibc' sideset on wp hole boundary.
- 'kelvin' material (+ Periodic identification) for open boundary.
- Canonical .jou: `src/radia/panels/samples/ih_fem_kelvin_skin_fine.jou`.

### Application interface

`IHDesignSpec` and the Induction Heating Simulink block present the above through
settings sections (Method / Drive / Coil material / Coil geometry /
Workpiece material / Workpiece impedance / Linear solver / Advanced).

Non-obvious features:

1. **Material presets** (Copper / Aluminum / Brass / Steel mu_r=100 /
   Steel mu_r=500 / Stainless 304 / Custom).  Selecting a preset
   auto-fills sigma+mu_r and disables their lineedits; `Custom`
   re-enables manual entry.  **Coil and wp materials are INDEPENDENT**
   — setting one does not change the other.

2. **Impedance model** (Linear SIBC / Nonlinear ESIM).  Selecting
   ESIM reveals bh_file, max_iter, tol widgets.  `calc_inductance`,
   `calc_fem_kelvin` and `calc_fem_coilmesh` run the ESIM Karl iteration
   for `--impedance-model esim`; it needs a 2-column BH file
   (`--bh-file`).  An unconverged ESIM cell solve always raises.  An
   unconverged outer Karl loop raises in `calc_fem_kelvin`, while
   `calc_inductance` and `calc_fem_coilmesh` report it as
   `"esim_converged": false` in their JSON — check that field.

3. **Linear solver** per-method:
   - PEEC+BEM: Dense LU (small) / HACApK (large, O(N log N))
   - HCurl FEM-SIBC: prefer AMS for p=1, BDDC with AMS coarse solve for
     p=2 or p=3 (`calc_fem_kelvin --solver auto`). The edge-only wirebasket
     uses three AMS coarse cycles per BDDC application. Both use COCR. Every
     route, direct included, has a true relative residual gate of 1e-6 on
     the original assembled free rows, with at most three residual corrections
     using the same operator. Iterations and correction residuals are recorded.
     ICCG is refused unless it converges.
   - Compound A-V: select `sparsecholesky` explicitly; the HCurl AMS
     validation does not establish support for compound or periodic spaces.
   - Scalar thermal H1: `sparsecholesky`; Maxwell AMS is not its solver.
   Direct FE solves use `sparsecholesky`, with no alternative direct fallback.
   The small P1-P3 SIBC parity lane validates accuracy, not a speed advantage.
   AMS and BDDC+AMS are not validated on periodic Kelvin: `auto` selects
   `sparsecholesky` there; explicit `ams` or `bddc` requests fail. Kelvin-interior
   nodal gradient modes are fixed by a boundary-rooted edge forest without
   changing curls or interface values. The residual gate includes the original
   rows removed from the factorization by that gauge.
   Filament loads are integrated separately between element transitions:
   integrating a closed loop across discontinuous basis gradients with one
   fixed quadrature rule can inject an unphysical gradient load.

4. **.vol label validation**: inspect `mesh.GetMaterials()` and
   `mesh.GetBoundaries()` against the selected solver's required contracts.
   Use the headless CLI validation errors; there is no `inspect_vol_labels`
   utility or retired analysis-panel Run-button contract to invoke.

5. **Physics sanity** shown under method:
   - wp delta and R/delta ratio; warns when R/delta < 3
     (SIBC approximation breaks, volumetric FEM preferred)
   - Post-solve: coil h_max vs delta (flags under-resolved coil mesh)

6. **Output** formatted in IH Summary section:
   - `L_coil (vacuum)` vs `L_total (with wp)` distinguished
   - `P_workpiece` (primary) / `P_coil (mesh-sensitive ±15%)` /
     `P_total` / `Heating efficiency = P_wp/P_total`
   - `H_t_rms`, `wp_area`, `coil delta`, `coil h_max` with OK/WARN tags

### Implementation notes (2026-04-19)

1. phi H1 DEFINED-ON coil only (`definedon=mesh.Materials('coil')`).
2. `A + grad(phi)` compound term in the eddy bilinear form on coil.
3. Dirichlet lift: set gf_phi = 1 on source, then solve
   `A x^(-1) (- A . gfu) + gfu = gfu_new`.
4. Post-solve volume current extraction uses a SEPARATE H1 scalar
   `psi_n` with dirichlet=sink, Set(1) on source, then
   `I_out = Integrate(J_coil * grad(psi_n) * dx('coil'))`.  Scale
   the entire gfu by I_target / I_out to normalize.
5. WP H_t from SIBC Robin:
   `H_t_rms = |jw/Z_s_wp| * sqrt(int|A_t|^2 dS / A_wp)`.
6. P_coil volumetric: `0.5 / sigma * int |J|^2 dV` over coil material.
   Warn if coil h_max > delta — mesh-sensitive accuracy.
"""

IH_PEEC_BEM_SIBC += SOLVER_CAPACITY


def get_ih_peec_documentation(topic: str = "all") -> str:
    """Return IH PEEC route documentation by topic."""
    topics = {
        "routes": IH_PEEC_ROUTES,
        "bem_sibc": IH_PEEC_BEM_SIBC,
    }
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join(topics.values())
    if topic in topics:
        return topics[topic]
    return f"Unknown topic: '{topic}'. Available: all, {', '.join(topics)}"
