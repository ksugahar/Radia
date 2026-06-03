# TEAM 28 levitation force via 3D (axisymmetric) Cauer Ladder Network

**Result (verified 2026-06-04):** a 6-stage Cauer Ladder Network (CLN)
reduced model reproduces the **TEAM Problem 28 electrodynamic-levitation
Lorentz force vs height** to better than 0.1%, and finds the levitation
equilibrium at **dZ = +4.1 mm** (full-FEM ~ +4 mm) where lift equals the
disk weight (~1.055 N).  To our knowledge this is the first time the CLN
reduction has been carried through to the actual TEAM 28 **levitation
force** (prior CLN-on-TEAM-28 work extracted only the decay spectrum
R_n / L_n on a generic test disk under a uniform field).

## What it shows

The coil-driven eddy-current disk is, at angular frequency `s = j*omega`,
a linear system `(K + s*N) X = F` where `K` is the s-independent
magnetostatic mixed phi-B operator, `N` the conductivity term
`v*(sigma*u/r)`, and `F` the coil source `v*Jz`.  The CLN / Cauer
reduction is the Krylov subspace generated from the **coil source** by the
magnetostatic-solve / sigma-accumulate recursion

    V_0     = K^{-1} F                          (s=0 coil response, no eddy)
    V_{k+1} = orthonormalise( K^{-1} (N V_k) )

The N-stage reduced model `(V^T K V + s V^T N V) y = V^T F`, evaluated at
50 Hz, gives the reduced field; the levitation Lorentz force follows from
the lab force integral `Fz = integral (Re[B_r]Re[J] - Im[B_r]Im[J]) 2*pi*r`
over the disk.  The force converges to the full-FEM value in ~5 stages:

| stages | F_z [N] | rel. err vs full |
|---|---|---|
| 1 | -0.047 | 97.8 % (DC, no eddy) |
| 2 | -2.278 | 3.9 % |
| 3 | -2.196 | 0.14 % |
| 5 | -2.19253 | 0.000 % |

## Files

| File | What |
|---|---|
| `team28_axisym_fem.py` | Repo-clean port of the lab full-FEM axisymmetric TEAM 28 solve (mixed phi-B + anisotropic-nu infinite shell). Reproduces the lab `.mat` force to **0.01%** at dZ=0. The ground-truth baseline. |
| `team28_cln_force.py`  | CLN/Cauer reduction at one height: builds K, N, F, shows the N-stage CLN force converging to full-FEM (golden). |
| `team28_cln_sweep.py`  | CLN force **vs height**, compared to the lab full-FEM `Fz1(dZ)`; recovers the levitation equilibrium ~+4 mm. |

## Source / provenance

- Geometry + full-FEM ground truth: lab learning material
  `W:\00_CAE\NGSolve\01_菅原\2024_08_TEAM28` (axisymmetric NGSolve TEAM 28,
  DC / 50Hz / 50Hz_可動 / Transient + field-validation figure).
  Disk: Al, R=65mm, t=3mm, sigma=3.4e7; coils: 960t/+20A (r=41mm) and
  576t/-20A (r=87.5mm) counter-wound, 50 Hz.
- CLN theory: `radia_mcp.mor` (mor_cln); Kameari-Ebrahimi-Sugahara-
  Shindo-Matsuo 2018, IEEE TMag 54(3):7201804.
- Method context: `radia_mcp.maglev` topics `cln_mor_control` /
  `radia_iem_fem`; the CLAUDE.md policy "Maglev Analysis: Radia + NGSolve".

## Run

```bash
python team28_axisym_fem.py   # full-FEM baseline  -> -2.1925 N @ dZ=0
python team28_cln_force.py    # CLN convergence    -> 5-stage golden
python team28_cln_sweep.py    # CLN force vs height -> equilibrium +4.1mm
```
