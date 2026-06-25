# axifem — NMR axisymmetric stored-reference validation

Reproduce a permanent-magnet NMR axisymmetric problem with the
Henrotte / Meeker P1-triangle pure-Python prototype, sample
`B_z(r)` along the symmetry plane, and **compare against a stored
reference and an NGSolve mixed-formulation reference at p=2** to
confirm the formulation tracks accepted axisymmetric solvers.

## Geometry

Cylindrical NMR magnet stack, axisymmetric:

- Domain: rectangle r ∈ [0, 70 mm], z ∈ [-50 mm, +50 mm]
- Magnet 1: r ∈ [0, 40 mm], z ∈ [-20, -10 mm], M = +z, H_c = 795 774.7 A/m
- Magnet 2: r ∈ [0, 40 mm], z ∈ [+10, +20 mm], M = +z, H_c = 795 774.7 A/m
- Air everywhere else
- Dirichlet φ = 0 on top, bottom, right, and on the axis

Sample `B_z` and `r·A_phi` at 1000 evenly-spaced points along
z = 0, r ∈ [0, 70 mm].

## How to run

The example bundles its own copy of `axifem_core.py` (the pure-Python
Henrotte/Meeker prototype) plus the two reference `.mat` files:

- `FEMM.mat` — solution from FEMM (David Meeker), p=2 reference
- `axisymmetric_mixed.mat` — solution from an NGSolve mixed-form
  axisymmetric solver, p=2 reference

so the script is fully self-contained:

```bash
cd examples/axifem/nmr_validation
python nmr_validate.py
```

Expected runtime: ~30–60 s on a typical workstation.  Outputs:

- `nmr_axifem_p1a.mat` — axifem B_z, B_r, A_phi, r·A_phi sampled
- `nmr_axifem_p1a.png` — three-method overlay plot:
  FEMM (black), NGSolve mixed (red dashed), axifem Phase 1a (blue)

## Expected output

The numerical comparison at sample points should look like

```
    r (mm)       B_z FEMM      B_z mixed      B_z axifem      ratio
      0.50   1.6028e-01   1.6029e-01     1.6014e-01    1.001
      2.00   1.6042e-01   1.6044e-01     1.6042e-01    1.000
      5.00   1.6114e-01   1.6114e-01     1.6107e-01    1.000
     10.00   1.6365e-01   1.6364e-01     1.6323e-01    1.003
     20.00   1.7028e-01   1.7028e-01     1.7020e-01    1.000
     40.00   3.660e-02    3.640e-02      5.008e-02     0.731
     60.00  -8.274e-02   -8.274e-02     -8.277e-02     1.000
```

`B_z stored ≈ B_z mixed ≈ B_z axifem` within ~0.3 % across r ∈ [0, 70 mm]
**except at r = 40 mm**, where all three solvers see the magnet edge
discontinuity in M_z and the P1-triangle axifem result diverges by
~27 %.  This is expected — high-order interpolation (Q2 quad) closes
that gap, see the disk_convergence example for the same lesson.

## Why this matters for IH

Axisymmetric magnets and induction-heating coils share the same
formulation: A_phi (or its FEMM-style proxy `phi = 2π r A_phi`) on
an (r, z) mesh, with material-region jumps in µ and σ.  The NMR case
exercises the **permanent-magnet RHS** branch (M_z dipole) of the
axifem assembly; the disk-convergence example exercises the
**eddy-current sigma-mass** branch.  Together they cover the two
operators that production IH simulations need.

The C++ port lives at `radia.axifem.AxiHenrotteFESpace`.
This pure-Python prototype is its **research-grade reference** —
when changes are made to the C++ side, the cross-validation test at
`tests/axifem/test_python_reference_consistency.py` ensures the
ported assembly still matches the prototype at machine precision.

## Source

Originally a research script under
`W:\30_CauerLadderNetwork\2026_04_01_長方形CLN\axifem`; the
reference `.mat` files originally lived under
`W:\00_CAE\NGSolve\01_菅原\2024_08_01_軸対称\NMR\`.  Both promoted
into this repository on 2026-05-10 with the axifem in-tree
integration so the example runs anywhere `pip install radia` (or a
clone) reaches.
