# radia-levitation

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
| `radia_levitation.mixed_galerkin.alpha` | Bulk Foster spectrum (NGSolve eigsh) + W(alpha) edge + Schur composition -> alpha(s) for any .vol |
| `radia_levitation.mixed_galerkin.cad_edges` | IGA-style CAD-direct edge extraction (12 edges exact for cuboid, mesh-independent) |
| `radia_levitation.ecb.lorentz` | F = integral J x B dV via Foster reconstruction, drag + lift on PM/coil source |
| `radia_levitation.ecb.plate_response` | alpha(s) frequency sweep + drag/lift crossover identification for plate ECB design |
| `radia_levitation.simulink.export` | State-space (A, B, C, D) export to MATLAB .mat + helper .m for Simulink LTI block |

## Quick start

```python
from radia_levitation.mixed_galerkin import alpha
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

| Example | Notebook / script |
|---------|-------------------|
| Single cube alpha(s) sweep | `examples/cube_alpha_sweep.py` |
| Linear ECB drag + lift vs velocity | `examples/ecb_linear.py` |
| Plate-over-PM levitation force | `examples/plate_levitation.py` |
| Simulink LTI export | `examples/simulink_export.py` |
| L-shape / multi-step bar | `examples/L_shape_response.py` |

## Force computation: Lorentz vs Maxwell stress

In this framework Lorentz force `F = integral J x B_ext dV` is MORE
ACCURATE than Maxwell stress tensor surface integral:

| Method | Error chain |
|--------|-------------|
| **Lorentz** | Foster truncation -> J(r), Radia analytical -> B_ext (exact), volume integral (1 step) |
| Maxwell stress | Foster -> J -> Biot-Savart -> B_induced, Radia -> B_ext, surface integral (3 steps) |

Reasons:
1. Maxwell needs an additional Biot-Savart from J (more numerical error).
2. Volume integration averages truncation error; surface integration
   does not.
3. Scalar Dirichlet formulation has known boundary artifacts that
   surface evaluation magnifies.

So the package defaults to Lorentz.  Maxwell stress is available
(`radia_levitation.ecb.maxwell_stress`) for cross-validation.

## Simulink integration

Export the polarizability as a continuous-time LTI state-space model:

```python
from radia_levitation.simulink import export
A, B, C, D = export.build_state_space(...)
export.save_mat("plant_lti.mat", A=A, B=B, C=C, D=D, ...)
```

Then in MATLAB:

```matlab
load('plant_lti.mat');
sys = ss(A, B, C, D);
bodeplot(sys);
% Simulink: place "LTI System" block -> set sys reference -> connect
%           to PM mechanics + controller in feedback loop
```

The state count is `n_foster + n_warburg + 1` (Foster bulk modes,
diffusive-quadrature Warburg rungs for K_SIBC/sqrt(s), one integrator
for c_1/s).  Typical state count: 50-300 depending on accuracy target.

## Related

- **radia-mcp** maglev knowledge: see `radia_iem_fem` and
  `cln_mor_control` topics (the published research line: Yano bachelor
  + Yano master, IEEE TMag 2018).
- **TEAM 28** electrodynamic levitation: `examples/team28/` reproduces
  the lab full-FEM force benchmark with 5-stage CLN reduction (0.000%
  error at the equilibrium height).
- **Sphere / ellipsoid / coil**: closed-form anchors and Radia open-
  boundary coil examples (migrated from
  `examples/CLN/scripts/team28_levitation/`).
