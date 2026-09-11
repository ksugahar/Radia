# radia.maglev

`radia.maglev` is the magnetic-levitation / eddy-current-brake
application module of **radia** (an application domain like `radia.ih`
(induction heating) or the electromagnet panel -- NOT a separately
released PyPI package; it ships inside the `radia` wheel and the
knowledge lives in the `radia_mcp.maglev` MCP subpackage).

Mixed Galerkin CLN-SIBC framework for magnetic levitation and eddy-current
brake analysis.  Given an arbitrary conductor geometry (`.vol` mesh or
OCC primitive), computes the polarizability tensor alpha(s) over a wide
frequency band, the Lorentz force F = integral J x B dV at any drive
condition, and exports a Simulink/MATLAB state-space LTI model for
control-coupled simulation.

## Theoretical basis

Mixed Galerkin reduction (Sugahara-Nagamine-Hane 2026, IEEE TMag in
prep): bulk CLN Krylov modes (s=0) Schur-composed with a polyhedral
surface envelope.  The polyhedral edge correction uses the closed-form
wedge function

    W(alpha) = (4/pi) cot(alpha/2)

(Wang-Lavers-Zhou 1992 surface impedance, pinned to the cuboid Mellin
anchor W(pi/2) = 4/pi).  See `journal_manuscript_paper1/master.tex` SIII
for the derivation.

## Modules

| Module | Purpose |
|--------|---------|
| `radia.maglev.mixed_galerkin.alpha` | Bulk Foster spectrum (NGSolve eigsh) + W(alpha) edge + Schur composition -> alpha(s) for any .vol |
| `radia.maglev.mixed_galerkin.cad_edges` | IGA-style CAD-direct edge extraction (12 edges exact for cuboid, mesh-independent) |
| `radia.maglev.ecb.lorentz` | Cycle-averaged 3-D HCurl-VIM plate force/torque through the common `radia.force` integrator; the older scalar Foster reaction-field approximation remains available with an explicit warning and direct-solve truncation fallback |
| `radia.maglev.ecb.plate_response` | alpha(s) frequency sweep + drag/lift crossover identification for plate ECB design |
| `radia.maglev.simulink.export` | State-space (A, B, C, D) export to MATLAB .mat + helper .m for Simulink LTI block |

## Quick start

```python
from radia.maglev.mixed_galerkin import alpha
from ngsolve import Mesh

mesh = Mesh("conductor.vol")
lam, tau, g_n, V = alpha.bulk_foster_via_eigen(mesh, sigma=5.8e7, mu=4e-7*3.14159, n_eigen=300)
S, edges = alpha.measure_total_area_and_edges(mesh)
K_SIBC = alpha.K_SIBC_total(S, 5.8e7, 4e-7*3.14159)
c1 = alpha.c1_polyhedral(edges, 4e-7*3.14159)

# alpha(s) at any frequency
import math
s = 1j * 2 * math.pi * 1e3   # 1 kHz
Y = alpha.Y_mixed(s, lam, tau, g_n, K_SIBC, c1)
a = alpha.alpha_from_Y(Y, V, sigma=5.8e7)
print(f"alpha(1 kHz) / V = {a/V:.4f}")
```

## Examples

Curated, runnable demonstrations live under `docs/maglev/demos/` (the broad
research corpus is in `validation_test/maglev/research_cln/`, see below):

| Example | Script |
|---------|--------|
| Single cube alpha(s) sweep (CAD-direct + mesh-derived) | `docs/maglev/demos/cube_alpha_sweep.py` |
| Sphere induced-dipole levitation force vs frequency | `docs/maglev/demos/sphere/maglev_sphere_force.py` |
| Coil-driven sphere eddy force | `docs/maglev/demos/sphere/coil_sphere_eddy_force.py` |
| Coil + sphere equilibrium height | `docs/maglev/demos/sphere/coil_maglev_equilibrium.py` |
| Ellipsoid demag / alpha tensor (axisym + 3D HCurl) | `docs/maglev/demos/ellipsoid/ellipsoid_alpha_tensor.py` |
| 3D cuboid CLN-SIBC standalone demo | `docs/maglev/demos/cuboid/cln_sibc_cuboid_3d.py` |
| TEAM 28 electrodynamic levitation (CLN force / sweep) | `docs/maglev/demos/team28/team28_cln_force.py` |

## Research corpus (`validation_test/maglev/research_cln/`)

`validation_test/maglev/research_cln/` holds the absorbed radia-cln research corpus — the raw
verification sweeps, Mathematica derivations, and iteration history that
back the theory but are **not** README-backed package examples (per the
lab Sample Promotion Ladder: tests -> examples -> panels).  Contents:

| Subdir | What |
|--------|------|
| `validation_test/maglev/research_cln/ngsolve_validation/` | NGSolve FEM verification sweeps (dated `*_results.json` snapshots, one-off probe scripts) |
| `validation_test/maglev/research_cln/axifem/` | Historical Henrotte-basis axifem prototype corpus; canonical artifacts now live under `docs/axifem/`, `tests/axifem/`, and `validation_test/axifem/` |
| `validation_test/maglev/research_cln/multiconn_loop_method/` | T-Omega multiply-connected bath-plate notebook |
| `validation_test/maglev/research_cln/tanimoto_canonical/` | Tanimoto canonical CLN notebooks (A-phi / A-T / T-Omega) |
| `validation_test/maglev/research_cln/*.wls` | Schur-F / CF / quadrupole / polarizability symbolic derivations |
| `validation_test/maglev/research_cln/bem_cln_*` | BEM-CLN multi-conductor verification (iteration history) |

Promote a script from `research/` to `docs/maglev/demos/` only after it gains a
README, runs standalone, and demonstrates one clear concept.

## Force computation: three-dimensional HCurl-VIM

The quantitative plate route is the sampled volume Lorentz integral
`F = integral J x B_ext dV`, with `J` supplied by the divergence-free 3-D
HCurl current basis and the open-boundary VIM interaction.  Use
`compute_lorentz_force_torque_via_hcurl_vim` for raw force/torque or
`compute_lorentz_force_result_via_hcurl_vim` for the common
`radia.force-result/v1` conductor/source pair.

The older scalar Dirichlet/Foster route is retained as a stated reduced model,
not as the physical reference.  Its direct scalar solve is useful for measuring
Foster truncation: `compute_lorentz_force_via_foster_verified` falls back to
that direct solve when a high-frequency basis is too short.  This does not cure
the scalar ansatz.  The independent rank- and mesh-converged 3-D lane
`validation_test/maglev/ecb_foster_lorentz_3d_reference.py` finds 63--89%
lift differences on its plate case and rejects the scalar model for
quantitative 3-D force.

There is not yet an ECB-specific enclosing-air Maxwell-stress route.  When one
is added, it should be treated as an independent extraction cross-check of the
3-D current solution, not as a correction applied to the scalar model.

## Simulink integration

The packaged `matlab/radia_maglev.slx` model is the human-facing dynamic
interface. Open it with `radia.simulink.openMagLev()`. Its masked MagLev plant
accepts `-dI/dt`, mechanical height, and coil current, then advances a
common-basis HCurl/CLN family and returns the induced port response and the
three-component Lorentz force. The tracked model carries diagnostic smoke data;
replace that data with a family exported by
`radia.vim.ExportHCurlEddyCLNFamilyJSON` for engineering use. Python is not
called during a Simulink time step.

For fixed-position diagnostics, export the polarizability as a continuous-time
LTI state-space model:

```python
from radia.maglev.simulink import export
A, B, C, D = export.build_state_space(...)
export.save_mat("plant_lti.mat", A=A, B=B, C=C, D=D, ...)
```

Then in MATLAB:

```matlab
load('plant_lti.mat');
sys = ss(A, B, C, D);
bodeplot(sys);
% Lower-level fixed-position plant for custom controller studies.
```

The state count is `n_foster + n_warburg + 1` (Foster bulk modes,
diffusive-quadrature Warburg rungs for K_SIBC/sqrt(s), one integrator
for c_1/s).  Typical state count: 50-300 depending on accuracy target.

## Related

- **radia-mcp** maglev knowledge: see `radia_iem_fem` and
  `cln_mor_control` topics (the published research line: Yano bachelor
  + Yano master, IEEE TMag 2018).
- **TEAM 28** electrodynamic levitation: `docs/maglev/demos/team28/` reproduces
  the lab full-FEM force benchmark with 5-stage CLN reduction (0.000%
  error at the equilibrium height).
- **Sphere / ellipsoid / coil**: closed-form anchors and Radia open-
  boundary coil examples (migrated from
  `docs/maglev/demos/{sphere,ellipsoid,team28}/`).
