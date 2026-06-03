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
| `cln_sibc_cuboid_3d.py` | Python port of the lab CLN-SIBC (Warburg-Schur) 3D cuboid core: Foster admittance + CLN reduction + Schur-F SIBC termination + polarizability `alpha(s)=V-Y/sigma`. The non-axisym building block. |
| `levitation_sphere_force.py` | **Isotropic induced-dipole AC levitation force** on a conducting sphere, coefficient pinned by the analytic perfect-conductor limit, frequency response reduced by CLN/Cauer. See below. |

## Isotropic levitation force (sphere) -- coefficient pinned, CLN-reduced

`levitation_sphere_force.py` builds the AC levitation force from the
induced magnetic dipole, with every constant verified.  A **sphere is
isotropic**, so the scalar polarizability already ported (`cln_sibc_
cuboid_3d.py`) applies directly -- no anisotropic tensor is needed to
demonstrate (and verify) a real levitation force.

The conducting sphere (Landau-Lifshitz ECM sec. 59) has magnetic response

    G(x) = -1/2 [ 1 - 3/x^2 + (3/x) cot x ],   x = (1+i) a / delta,

with `G(0)=0` (DC, no eddy response) and `G(inf)=-1/2` (perfect-conductor
flux exclusion).  The time-averaged levitation force on the induced
dipole in a field gradient is

    <F> = (pi a^3 / mu0) Re[G(x)] grad(B0^2),     Re[G] < 0  =>  LIFT.

| check | result |
|---|---|
| limits of G | DC `Re G -> 0`, HF `Re G -> -0.4997` |
| sign | `Re G < 0` for all f in [1 Hz, 100 MHz] -> lift at every frequency |
| CLN/Cauer reduction | stage 4 within **0.013%**, stage 6 **0.0000%** of the full modal system |
| coefficient pin | HF lift `31.22 mN` vs perfect-conductor `(pi a^3/2 mu0)|grad B0^2| = 31.25 mN` (**0.09%**) |

The same `(pi a^3 / 2 mu0) grad(B0^2)` coefficient is derived independently
from the perfect-conductor energy `U = -1/2 m.B` and reproduced by the
induced-dipole formula -- so the complex-AC sign and normalization are
pinned, not guessed.  The lift rises from ~0 (DC) through the eddy-current
transition (`a/delta ~ 1-5`) to the perfect-conductor saturation, exactly
as expected.  Isotropic; the cuboid `a!=b!=c` alpha tensor is a separable
refinement (not required for the force).

**Note -- "anisotropy" here = SHAPE, not material.**  The cuboid
`alpha = diag(alpha_x, alpha_y, alpha_z)` is direction-dependent because
the dimensions `a, b, c` differ, not because the conductor is an
anisotropic material -- copper stays a scalar `sigma`/`mu`.  A field along
`z` drives eddy currents in the `a x b` cross-section, along `x` in the
`b x c` cross-section, etc., so `a != b != c` gives three different eddy
time constants and hence a direction-split response.  It is the AC /
eddy-current generalization of the magnetostatic **demagnetizing-factor
tensor** (sphere: isotropic `1/3`; ellipsoid / brick: direction-dependent,
from shape alone).  Material anisotropy (tensor `sigma` / `mu`) is a
separate, genuinely-material effect, not what this refinement is about.

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
python team28_axisym_fem.py      # full-FEM baseline  -> -2.1925 N @ dZ=0
python team28_cln_force.py       # CLN convergence    -> 5-stage golden
python team28_cln_sweep.py       # CLN force vs height -> equilibrium +4.1mm
python cln_sibc_cuboid_3d.py     # CLN-SIBC 3D cuboid core (alpha, Schur-F)
python levitation_sphere_force.py  # isotropic levitation force, coeff pinned
```
