# User guide — Equivalence-theorem near-field source

This is the user-facing guide for `radia.equivalence_source.NearFieldSource`.
For the internal C++ architecture see [`CPP_DESIGN.md`](CPP_DESIGN.md).
For Phase D acceleration plans see [`FMM_DESIGN.md`](FMM_DESIGN.md).

## 1. What it does (one paragraph)

Given the magnetic field **H** (and optionally electric field **E**)
on a **closed surface** enclosing some source region, reconstruct
the field at any external observation point via the
**Stratton-Chu surface integral** (Schelkunoff / Love equivalence
theorem).  Use it when you want to:

- Save a FEM-solved interior field as a portable file and re-evaluate
  the **exterior** field at arbitrary points without redoing the FEM
  solve.
- Couple a Radia / NGSolve interior solve to an external MoM /
  ray-tracing tool via the Nastran NFS format.
- Verify the **equivalence-theorem null-field property**: evaluation
  inside the surface returns ≈ 0, exterior returns the true field
  (see `validation_test/equivalence_source/null_field_property.py`).

The implementation lives at
[`src/radia/equivalence_source.py`](../../src/radia/equivalence_source.py)
with a C++ accelerator at
[`src/core/rad_equivalence_source.cpp`](../../src/core/rad_equivalence_source.cpp).

## 2. When to use which kernel

| Regime | Kernel | Method |
|---|---|---|
| **Magnetostatic** (ω = 0): permanent magnets, DC coils, soft iron | `1/R` Laplace | `evaluate_static_H(obs)` — Phase A C++ |
| **MQS / Darwin** (low ω, kHz): induction heating, eddy current | dyadic full-wave | `evaluate(obs, omega)` — Phase B C++ |
| **Time-harmonic full-wave** (MHz–GHz): WPT, EMC, antennas | dyadic full-wave | same `evaluate(obs, omega)` |
| **Full radiation** (far field) | classical Stratton-Chu | same `evaluate(obs, omega)` |

The static and harmonic C++ kernels are both verified production paths.
The harmonic dyadic path is tracked by
`validation_test/equivalence_source/phase2_wpt_harmonic.py`; the
2026-06-28 LAB run reconstructs the 1 MHz Hertzian-dipole case within
0.12% for both E and nonzero H observations.

## 3. Quickstart

```python
import numpy as np
from radia.equivalence_source import NearFieldSource

# 1. Build a closed surface (here: sphere triangulation around the origin)
centroids, normals, areas = NearFieldSource.spherical_surface_mesh(
    R=0.5, n_theta=30, n_phi=60)

# 2. Sample H on that surface (here: analytical dipole, mz = 1 A·m²)
def H_dipole(p):
    r = np.linalg.norm(p); rhat = p / r
    return 1.0 / (4*np.pi*r**3) * (3 * rhat[2] * rhat - np.array([0., 0., 1.]))

H_surface = np.array([H_dipole(c) for c in centroids])

# 3. Construct NearFieldSource (E=None for pure magnetostatic)
nfs = NearFieldSource.from_surface_mesh(
    centroids, normals, areas, E=None, H=H_surface)

# 4. Evaluate H at exterior obs points -- ONE-LINE, C++ accelerated
obs = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 2.0]])
H_reconstructed = nfs.evaluate_static_H(obs)
# H_reconstructed.shape == (3, 3): rows = obs, cols = (Hx, Hy, Hz)
```

For time-harmonic (e.g. 13.56 MHz WPT):

```python
omega = 2 * np.pi * 13.56e6
# Need BOTH E and H on the surface for harmonic
nfs = NearFieldSource.from_surface_mesh(centroids, normals, areas,
                                          E=E_surface, H=H_surface,
                                          omega=omega)
E_rec, H_rec = nfs.evaluate(obs, omega=omega)
# Both are complex (M, 3) arrays.
```

## 4. Constructors — pick the one that fits your source

### 4.1 From analytical or hand-rolled samples — `from_surface_mesh`

```python
nfs = NearFieldSource.from_surface_mesh(
    centroids,   # (N, 3) face centroids [m]
    normals,     # (N, 3) OUTWARD unit normals
    areas,       # (N,)   face areas [m²]
    E=E_complex_or_real,    # (N, 3) or None for static
    H=H_complex_or_real,    # (N, 3) -- required
    omega=0.0,              # angular frequency [rad/s]
)
```

### 4.2 From a NGSolve FEM solve — `extract_ngsolve`

```python
nfs = NearFieldSource.extract_ngsolve(
    mesh,           # NGSolve mesh (must have the surface as a BND label)
    gf_E=gf_E,      # NGSolve GridFunction(HCurl) on the volume — optional for static
    gf_H=gf_H,      # NGSolve GridFunction(HCurl) on the volume — required
    surface_label="nfs_surface",   # mesh boundary name to extract on
    omega=2*np.pi*60e3,
)
```

The boundary surface must be a CLOSED surface in the mesh, labelled
in the `.vol` (`SetBCName(...)` or via a Cubit sideset named
`nfs_surface`).  Sampling is centroid-based per face.

### 4.3 From a Radia container — `extract_radia`

```python
rad.UtiDelAll()
container = build_pm_array(...)   # rad.ObjCnt of permanent magnets
nfs = NearFieldSource.extract_radia(
    container,
    surface_skeleton=(centroids, normals, areas),   # e.g. from spherical_surface_mesh
)
# H is sampled from rad.Fld() at each centroid; static only.
```

### 4.4 From a `.nfs.json` file — `load`

Round-trip persistence:

```python
nfs.save("my_source.nfs.json")
# ... later, or on another machine ...
nfs = NearFieldSource.load("my_source.nfs.json")
H_ext = nfs.evaluate_static_H(obs)
```

The JSON format is human-readable, small (panel count × 9 floats
typically <1 MB for 10⁴ panels), and version-tagged.

## 5. Workflow: Cubit → .vol → NGSolve → NearFieldSource → external eval

This is the canonical production pipeline; see
[`validation_test/equivalence_source/phase3_e2e_cubit_to_sol.py`](../../validation_test/equivalence_source/phase3_e2e_cubit_to_sol.py)
for a verified end-to-end test.

```
┌─────────────────────────────────────────────────────────────────┐
│  Cubit headless: build mesh, label "nfs_surface" sideset        │
│  export netgen "model.vol" order 2 overwrite              │
└──────────────────────────┬──────────────────────────────────────┘
                           │ .vol with bcnames "nfs_surface"
┌──────────────────────────┼──────────────────────────────────────┐
│  NGSolve FEM solve                                              │
│    mesh = Mesh("model.vol")                                     │
│    # FEM-Kelvin / Omega-reduced / whatever                      │
│    gf_H = ...   (HCurl GridFunction)                            │
│    gf_H.Save("model_H.sol")                                     │
└──────────────────────────┬──────────────────────────────────────┘
                           │ gf_H on "nfs_surface" BND
┌──────────────────────────┼──────────────────────────────────────┐
│  nfs = NearFieldSource.extract_ngsolve(                         │
│             mesh, gf_H=gf_H,                                    │
│             surface_label="nfs_surface", omega=0)               │
│  nfs.save("model.nfs.json")                                     │
└──────────────────────────┬──────────────────────────────────────┘
                           │ portable .nfs.json
┌──────────────────────────┼──────────────────────────────────────┐
│  Anywhere -- different machine, different solver                │
│    nfs = NearFieldSource.load("model.nfs.json")                 │
│    H_far = nfs.evaluate_static_H(my_obs_grid)                   │
│  or: rad_objs = nfs.to_radia_objects()                          │
│       feed into another Radia simulation as a background field  │
└─────────────────────────────────────────────────────────────────┘
```

## 6. Performance

| Path | Speed | Use when |
|---|---|---|
| Phase A C++ `evaluate_static_H(use_cpp=True)` | typically tens of times faster than Python at N=10³–10⁵ | **Default** |
| Phase B C++ `evaluate(use_cpp=True)` | similar (40–80× harmonic) | ω > 0 default |
| Python fallback (`use_cpp=False`) | slow (10⁵-panel × 100-obs ≈ 30 s) | Regression cross-check ONLY |

Benchmarks (see `validation_test/equivalence_source/bench_static.py`):
- 2026-06-28 LAB run: 4/4 numerical cases passed.
- Speed target result: 4/4 cases met the production-scale 50× target.
- Bit-identical results: `‖H_C++ − H_python‖∞ ≈ 1e-15` (rounding noise).

**Above N_face × N_obs > 10⁹**: direct C++ becomes the bottleneck.
See [`FMM_DESIGN.md`](FMM_DESIGN.md) for the NGSolve.bem-based
Phase D acceleration plan.

## 7. Adapting accuracy

### 7.1 Panel density

The reconstruction error scales as **1/N** in panel count (for a
spherical extraction surface, this is the discretisation error of
the surface integral with a single centroid quadrature point per
face).  Doubling `n_theta × n_phi` halves the error.

Typical numbers (Phase 1 magnetic dipole golden):
- N = 14 160 panels → 0.83% max rel error
- N = 56 640 panels → 0.2% max rel error
- N = 100 000 panels → 0.1% max rel error

### 7.2 Extraction-surface radius

The surface only needs to **enclose all sources**.  Placing it too
close to a source increases the sensitivity to per-face sampling
error; placing it too far costs panels (since you typically want
roughly uniform angular resolution).  Rule of thumb: place at
1.5–3× the source diameter.

### 7.3 Static vs harmonic — when does ω = 0 stop being OK?

`evaluate_static_H` is the Phase A magnetostatic reduction.  It is
**exact for ω = 0** (DC) — no approximation beyond the discretisation
error of the surface integral.

For ω > 0 use `evaluate(obs, omega=...)`.  The dyadic-correct C++
kernel (Phase B) matches Phase A static at ω → 0 to 5e-15 relative
error (machine precision) — so the harmonic call is always safe at
low frequency, but the static call is faster (real-only arithmetic,
no exponentials).

Threshold (informal): if `ω · R_obs / c < 1e-6`, use static — it's
~2× faster.  Otherwise use harmonic.

## 8. Limitations

- **Single quadrature point per face** (centroid).  Higher-order
  Gauss quadrature is a roadmap item (estimated 2× accuracy at
  same panel count, ~1.5× cost).
- **Static reconstruction is H-only** — no E reconstruction in the
  Phase A path (would just give zero).  For E set `omega > 0` and
  use `evaluate(...)`.
- **Surface must be closed and enclose the source**.  The
  equivalence-theorem null-field property fails on open / partial
  surfaces.  See `validation_test/equivalence_source/null_field_property.py`
  for a verification.
- **Obs points must be outside the enclosing surface** for physical
  results.  Inside, you get ≈ 0 (the null-field property), NOT the
  interior field of the original source.  This is BY DESIGN; the
  interior field requires the original FEM mesh.

## 9. See also

- [`CPP_DESIGN.md`](CPP_DESIGN.md) — C++ kernel architecture
- [`FMM_DESIGN.md`](FMM_DESIGN.md) — Phase D acceleration via NGSolve.bem
- `validation_test/equivalence_source/` — executable validation corpus
  - `phase1_static_coil.py` — golden test: analytical coil
  - `phase2_wpt_harmonic.py` — harmonic Hertzian dipole, 1 MHz deep-near-field PASS
  - `phase3_e2e_cubit_to_sol.py` — Cubit → NGSolve → NFS full e2e
  - `null_field_property.py` — interior null + exterior reconstruction
  - `bench_static.py` — Phase A C++ benchmark
- MCP tool: `radia_mcp.fem.equivalence_source_knowledge`
- Sugahara Lab 2008 axi slides — original derivation
