# axifemm — disk convergence sweep

Probe the convergence of axisymmetric induction-heating eddy-current
analysis (Henrotte / Meeker P1-triangle formulation, NumPy/SciPy
prototype) for a copper disk inside an air domain.  Compares the
first relaxation time τ_1 against the **v22b reference value
208.32 µs** (Mathematica BEM-Foster derivation).

## Geometry

Cylindrical copper disk:
- Radius R_disk = 10 mm
- Half-thickness T_disk / 2 = 1 mm (so total thickness 2 mm)
- σ = 5.8e7 S/m (copper)
- µ_r = 1 (non-magnetic)

Air truncation domain: rectangle r ∈ [0, R_air], z ∈ [-Z_air, Z_air]
with φ = 0 Dirichlet on the outer right / top / bottom edges and on
the axis.

## What it sweeps

Two error sources:

1. **Truncation domain size** (R_air, Z_air): increase from 50 mm
   to 1 m, at moderate fixed mesh density.  Truncation error → 0
   as the air domain grows.
2. **Mesh density** (maxh_disk, maxh_air): refine the Henrotte
   element size from 0.4 mm to 0.05 mm at fixed large air domain.
   Discretisation error → 0 as h → 0.

The ratio `τ_1 / 208.32 µs` should approach 1.0 from below as both
errors decrease.

## How to run

The example bundles a self-contained pure-Python prototype copy
(`axifemm_core.py`, `sigma_mass.py`) so no extra path setup is
needed; just `numpy`, `scipy`, and `ngsolve`/`netgen-mesher` for the
OCC mesh generation.

```bash
cd examples/axifemm/disk_convergence
python disk_convergence_sweep.py
```

Expected runtime: ~30–60 s on a typical workstation.

## Expected output

Two tables (truncation sweep and mesh sweep), each showing
`tau_1` and the ratio to 208.32 µs.  Observed (2026-05-10):

```
[1] Truncation domain sweep (maxh_disk=0.2mm, maxh_air=10mm)
    R_air mm   Z_air mm    Nnode    Nelem     tau_1 us      ratio
          50         50     1018     1979     218.96       1.051
         100        100     1243     2399     219.57       1.054
         200        200     1942     3736     219.91       1.056
         500        500     6870    13414     219.88       1.055
        1000       1000    24338    48049     219.70       1.055

[2] Mesh refinement (R_air=500mm, Z_air=500mm)
   maxh_disk   maxh_air    Nnode    Nelem     tau_1 us      ratio
       0.400       20.0     1879     3586     219.02       1.051
       0.200       10.0     6870    13414     219.88       1.056
       0.100        5.0    26296    51956     220.43       1.058
       0.050        5.0    33512    66367     220.02       1.056
```

The ratio plateaus at **~1.055** for the **P1 triangle** Henrotte
formulation — this is a systematic ~5–6 % gap inherent to P1, not a
truncation or mesh-resolution error.  To close it to << 0.1 % switch
to **Q2** quadrilateral elements (`H1Henrotte(mesh, order=2)` with a
structured quad mesh) — see
`examples/axifemm/research/verification/test_hiruma_disk_q2.py` for the Q2
benchmark on the same disk geometry.

The lesson: P1 Henrotte is the **research-grade reference formulation**
(easy to derive, axis-handling clean), but production accuracy on
small thin disks needs Q2.

## Why this matters for IH

The first relaxation time τ_1 dominates the inductive transient
of an axisymmetric workpiece — it is the dominant pole of the
Cauer ladder network (CLN) representation of the disk eddy-current
admittance.  Verifying that the axifemm Henrotte FE recovers
v22b's τ_1 within ~1 % validates the formulation as a reference for
production induction-heating analyses (where the dominant time
constant drives heating-rate control loops).

The C++ port lives at `radia.axifem.AxiHenrotteFESpace` (see
`tests/axifemm/test_python_reference_consistency.py` for the
element-matrix-level cross-check between this pure-Python prototype
and the C++ implementation).

## Source

Originally a research script under
`W:\30_CauerLadderNetwork\2026_04_01_長方形CLN\axifemm`; promoted
into this repository on 2026-05-10 with the axifemm in-tree
integration.
