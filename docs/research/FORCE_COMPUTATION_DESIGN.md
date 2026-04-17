# Force / Torque / Energy Computation: Radia-NGSolve Design

**Status**: Research phase. No implementation planned yet. Thorough design required before coding starts.

## Background

The legacy `rad.FldEnr` / `rad.FldEnrFrc` / `rad.FldEnrTrq` API was removed in Phase C
(commit `98d7f65`, 2026-04-16) along with the `SubdivideItself*` family of methods.

The removed API used a 1990s-era approach: physically subdivide the destination magnet
into `kx * ky * kz` pieces (user-controlled via `SbdPar`) and apply a midpoint
quadrature rule to the pair-energy integral. Disadvantages:

- Physical mesh modification (memory growth, cache thrashing).
- Fixed integer subdivision; no adaptive error control.
- `SubdivideItself*` alone accounted for ~3000 lines in `rad_polyhedron.cpp`.
- Midpoint rule: low-order convergence; dense subdivision is required for engineering
  accuracy.

This document captures the design direction for a modern replacement.

## Core story

> **Using Radia + NGSolve, you eliminate the air mesh AND obtain stable
> high-accuracy electromagnetic force.**

| Aspect                | Pure FEM (NGSolve alone)           | Radia + NGSolve                        |
|-----------------------|-------------------------------------|----------------------------------------|
| Air mesh              | Required (large)                    | **Not required** (Biot-Savart + MSC)   |
| Open boundary         | PML or Kelvin transform             | **Natural open / Kelvin transform**    |
| Force calculation     | Nodal force (noisy) / Maxwell tensor (surface-dependent) / virtual work (two analyses) | **Radia analytical pair + NGSolve high-order quadrature** |
| Discretization error  | O(h^p) everywhere                   | **Zero on pair kernel; high-order Gauss on the rest** |
| Reproducibility       | Mesh-dependent                      | **Mesh-independent once converged**    |

## Why Radia is suited for this research

- **Built-in accuracy-verification environment**:
  `examples/ngsolve_integration/verify_curl_A_equals_B/` and related cross-validation
  tests already exist. Any new force kernel can be validated quantitatively in the
  same framework.
- **High-order element support via NGSolve integration**: Netgen 6.2.2603 provides
  hex/prism curving, `mesh.Curve(p)` for p = 1..5, and `Integrate(..., order=N)`
  with high-order Gauss quadrature. Radia's analytical kernel composed with NGSolve's
  high-order quadrature gives a unique position where discretization error can be
  controlled from both sides.
- **A combination no other software offers**: Pure FEM (COMSOL, Ansys, NGSolve alone)
  needs an air mesh, and force is always noisy (nodal / Maxwell / virtual work). Pure
  BEM (the legacy Radia `FldEnr`) was capped by `SubdivideItself` midpoint quadrature.
  The Radia-NGSolve hybrid is the only path that delivers "no air mesh × stable high
  accuracy" for force research.

## Integration methods, in order of preference

### 1. Radia analytical pair-interaction (reuse MSC/MMM kernel)

The existing kernels `Compute3x3BlockFast` (MMM tet-tet) and `Compute5x5/6x6BlockFast`
(MSC wedge / hex) return a closed-form demagnetization tensor `N_{ij}` between element
pairs. Derived quantities:

- Pair energy:       `W_{ij} = sigma_i . N_{ij} . sigma_j` (MSC)
                    or `m_i . N_{ij} . m_j` (MMM)
- Force on body i:   `F_i = sum_j (dN_{ij} / dr) . (sigma_i × sigma_j)`
                    -- analytical position-derivative of the kernel

**Zero discretization error**, no mesh subdivision, ~500 LOC (kernel reuse).

### 2. NGSolve `Integrate()` with high-order quadrature

Use case: Radia MSC/MMM solution is the source, an NGSolve mesh is the destination
(or vice-versa -- FEM domain coupled to Radia source).

```python
W = Integrate(M_dst * B_src, mesh_dst, order=5)  # 5-point Gauss on tets -> ~1e-8 accuracy
F = ...                                          # finite-diff on W, or Maxwell tensor on a boundary
```

`B_src = rad.RadiaField(src_obj, 'b')` is already provided as an NGSolve
`CoefficientFunction` (`_radia_pybind.pyd`). No new C++ is needed on that side.

### 3. Maxwell stress tensor (extended `rad.FldFrc`)

Today `rad.FldFrc` integrates the Maxwell tensor over a rectangular surface. Extension:
support OCC shapes or NGSolve boundaries. Still surface-placement dependent, but stable
when the surface is reasonably far from the magnet.

### 4. FEM nodal force / virtual work

**Not adopted.** Both are noisy and mesh-dependent.

## Proposed API (not implemented; requires review)

```python
# Option 1 -- Radia analytical pair-interaction (closed-form, zero error)
F = rad.AnalFrc(dst_obj, src_obj)
W = rad.AnalEn(dst_obj, src_obj)
T = rad.AnalTrq(dst_obj, src_obj, P)

# Option 2 -- NGSolve-based (when an NGSolve mesh is available)
from radia_ngsolve import force_from_energy, force_from_tensor
F = force_from_energy(dst_mesh, src_obj, order=5)                # volume integral
F = force_from_tensor(dst_mesh, dst_boundary_name, src_obj)      # boundary Maxwell tensor
```

The legacy `rad.FldEnr(dst, src, SbdPar)` is gone and will not be resurrected. The
`SbdPar` argument has no place in the new design.

## Design questions to resolve before implementation

1. **Position-derivative of the pair kernel**: extend `Compute*BlockFast` to
   `Compute*BlockFastGrad`. Closed-form derivatives of the solid-angle kernel exist,
   but the code generation is large. Consider a Mathematica / SymPy pass that emits
   `.cpp`.
2. **Body-type combinations**: PM × PM, PM × soft iron, soft iron × soft iron, PEEC ×
   magnetic. Per-combination API or one unified entry point?
3. **Torque reference point**: world origin, destination centroid, or user-specified?
   (Legacy API took `double* pP`.)
4. **NGSolve coupling recipe**: worked examples for `RadiaField` + `Integrate(order=N)`
   do not yet exist. Add to `examples/ngsolve_integration/`.
5. **Maxwell-tensor extension**: accept arbitrary closed surfaces (OCC shape -> surface
   quadrature) for `rad.FldFrc`.
6. **Test strategy**: analytical configurations (Smythe sphere, two PMs on an axis,
   coil-iron pair) with required relative error < 1e-6.
7. **Benchmark strategy**: compare against FEM force methods (nodal / Maxwell /
   virtual work) on mesh-refinement sweeps. Quantify "stability" as the noise envelope.

## Legacy code reference

The removed code is still reachable in git history. Use it as a reference *for the
math*, not as a template:

```bash
git show 98d7f65^:src/core/rad_geometry_3d.cpp  # ActualEnergyForceTorqueComp body
git show 98d7f65^:src/core/rad_geometry_3d.cpp  # EnergyForceTorqueCompAutoDestSubd
git show 98d7f65^:src/core/rad_geometry_3d.cpp  # ProceedNextStepEnergyForceTorqueComp
git show 98d7f65^:src/lib/radentry.cpp          # RadFldEnr / RadFldEnrFrc / RadFldEnrTrq
```

`ValidateForceChar` / `ValidateTorqueChar` parsed IDs of the form `"fx"`, `"fy"`,
`"fz"`, `"tx"`, `"ty"`, `"tz"`. The new API should use structured numpy-array returns
rather than ID strings.

## Status and next step

- Phase C complete (`98d7f65`, on `main`).
- Force API today: only `rad.FldFrc` (Maxwell tensor on a rectangular surface).
- **Before the next implementation session, re-read this document and resolve the 7
  design questions above.**
- User is not in a rush. Thorough design is the priority.
