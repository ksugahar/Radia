"""TEAM force/motion + nonlinear magnetostatic + eddy current problems: 2(cyl), 6, 7, 13, 17, 20, 21, 21a-1, 23, 28, 33b — lab core."""

OVERVIEW = r"""
# TEAM force / motion / nonlinear-B / eddy-current problems

| # | Name | Physics class | Lab use |
|---|------|--------------|---------|
| **2†** | Conducting Cylinder in Uniform Axial B ★ | A-Φ Dirichlet box, scattered field, I₀ Bessel ref | lab production |
| **6**  | Conducting Sphere Shell in Uniform B ★ | A-Φ Dirichlet box, scattered field, analytical ref | lab production |
| **6s** | Solid Conducting Sphere in Axial B ★ | A-Φ, sinh/cosh spherical Bessel ref | lab production |
| **6m** | Permeable Sphere (Magnetostatic) ★ | scattered-field A-form, demagnetization factor | lab production |
| **7**  | 3-D Eddy Current ★ | time-harmonic A-Φ, eddy current | lab production |
| **13** | 3-D Nonlinear C-Core B-Flux ★ | nonlinear magnetostatics, B-flux | lab production |
| 17 | Jumping Ring | Lorentz transient | demo / education |
| **20** | 3-D Static Force ★ | Maxwell stress, magnetostatic | lab production |
| **21a-0** | Tank-Wall Eddy Current ★ | A-Φ Dirichlet box, eddy loss | lab production |
| **21a-1** | Tank-Wall + 1 Central Slit ★ | same, slit blocks main current loop | lab production |
| **21b-2N** | Tank-Wall 2 Non-mag Plates + Gap ★ | same coils, 2×248mm plates with 200mm gap | lab production |
| **23** | Forces in Permanent Magnets ★ | PM force | lab production |
| 28 | Electrodynamic Levitation | levitation, force balance | maglev research |
| **33b** | Electric Local Force ★ | local force formulations | per-node verification |

★ = lab core validation problems.  † = TEAM-2 analog (exact spec from PDF, password protected).
"""


PROBLEM_7 = r"""
# TEAM Problem 7 — 3-D Eddy Current in Aluminium Plate (★ lab A-Φ formulation test)

## Geometry (Nakata/Takahashi/Yoshida, COMPUMAG 1988)

A thick aluminium plate with a rectangular hole, driven by a 50 Hz rectangular
coil. Eddy currents form in the plate; measurements: Bz on two horizontal lines
above the plate (A1-B1 at z=34mm, y=72/144mm) and Jy through the plate (A3/A4).

```
Plate : 294×294×19 mm, Al (σ=35.26 MS/m), hole 108×108 mm at offset (18,18) mm
Coil  : rectangular, 100×100 mm cross-section (wall 25 mm), 100 mm height,
        2742 turns, NI≈1000 A (coil centred at ~(394,100) mm in xy, z=49..149 mm)
Freq  : 50 Hz  (skin depth δ = 1/√(π·f·σ·μ₀) ≈ 12 mm in Al — partially penetrating)
```

## Open boundary: Kelvin transform (inner + outer sphere, Periodic BC)

Inner sphere: centre (147,147,50) mm, radius R=300 mm, box-clipped to remove 2 caps
  clip box: (cx-280, cy-280, cz-600) to (cx+600, cy+600, cz+600) mm
Outer sphere: same clipped shape, shifted +2R=600 mm in +x
Periodic BC: `ext_dome.Identify(int_dome, "periodic", IdentificationType.PERIODIC)`
GND vertex at outer sphere centre → Φ=0 gauge for the H1 space.
outer_cx = sphere_cx + 2*R = 0.147 + 0.60 = 0.747 m

## A-Φ formulation (weak, e^{jωt} convention)

State: A ∈ Periodic(HCurl(complex=True, nograds=True)),
       Φ ∈ Periodic(H1(complex=True, dirichlet_bbbnd="GND"))

Bilinear form (5 coupling + 1 regularisation):
  ∫ ν curl A·curl v dx
  + jωσ [A·v + ∇Φ·v + A·∇ψ + ∇Φ·∇ψ] dx_plate     ← all 4 A-Φ coupling blocks
  + ε ν ∇Φ·∇ψ dx                                     ← Φ null-space regularisation (ε=1e-10)

Outer-domain reluctivity: ν_outer = (r'/R)² ν₀   (r' = dist from outer sphere centre)

## Coil current source (tau_coil / pot_coil)

Inner coil square (projected): x ∈ [0.144, 0.244], y ∈ [0.050, 0.150]
  projx = clamp(x, 0.144, 0.244),  projy = clamp(y, 0.050, 0.150)
  tau_coil = (projy-y, x-projx, 0) / |tau|   [unit winding tangent]
  pot_coil = (0, 0, sqrt((projy-y)²+(projx-x)²) - 0.050)

Linear form (N=2742, L=0.100 m, dw=0.025 m):
  lf += -N/L · tau_coil · v_A.Trace() · ds("coili", bonus_intorder=4)
  lf += N/(dw·L) · pot_coil · curl(v_A) · dx("coil", bonus_intorder=4)

## NGSolve recipe (REUSE THIS — shipped as solve_eddy_current_harmonic_APhi)

```python
gfAPhi = solve_eddy_current_harmonic_APhi(
    mesh, nu, sigma, omega, add_source,
    order=5,         # order=3 is ~4x faster, ~20 % less accurate
    precond="local"  # "bddc" fails for complex=True (Cholesky NaN)
)
gfA, gfPhi = gfAPhi.components
B = curl(gfA)
J = -1j * omega * sigma * (gfA + grad(gfPhi))   # eddy current density [A/m²]
```

GMRes params (validated): tol=1e-8, maxsteps=2000, restart=200.
  order=3: ~180 iters / ~22 s.   order=5: ~200 iters / ~100 s.

## Reference values (TEAM-7 measured, Gauss, complex)

Line A1-B1 (y=72mm, z=34mm, x=0,18,...,288 mm — 17 points):
  Real: -4.9, -17.9, -22.1, -20.2, -15.7, +0.4, +43.6, +78.1, +71.6, +60.4,
        +53.9, +52.6, +53.8, +56.9, +59.2, +52.8, +27.6  [Gauss]
  Imag: -1.2, +2.5, +4.2, +4.0, +3.1, +2.3, +1.9, +5.0, +12.6, +14.2,
        +13.0, +12.4, +12.1, +12.3, +12.7, +10.0, +2.3  [Gauss]

Line A2-B2 (y=144mm, z=34mm): peak real ≈ +60 G at x≈144mm.

NGSolve <-> reference sign convention:
  Bz_ref_real = -1e4 * B[2].real [T→G conversion + sign flip]
  Bz_ref_imag = +1e4 * B[2].imag

## Cross-validation status (radia-mcp, 2026-06-04)

``solve_eddy_current_harmonic_APhi`` implemented and shipped.
Regression test: ``test_team7_eddy_bz`` (order=3, mesh_maxh=0.05 m, RMSE<12 G on A1-B1).
Yano's production validation: order=5, same mesh, visual match in comparison_BzJy.png.

## Lab source

Reference Python solver + output plots:
  W:\\00_CAE\\NGSolve\\矢野\\2026_03_30_TEAM_benchmark\\Problem7\\APhi\\
  Problem_7_A_Phi_Kelvin.py
  output_A_Phi_Kelvin_FreqDomain\\comparison_BzJy.png
  output_A_Phi_Kelvin_FreqDomain\\Bz_comparison.csv
"""


PROBLEM_20 = r"""
# TEAM Problem 20 — 3-D Static Force Problem (★ lab core)

## Geometry (official — Nakata/Takahashi/Morishige, Okayama Univ.)

A vertical STEEL center pole + a STEEL yoke (return frame), a DC coil of 381
turns wound around the pole, and an air gap between the pole top and the yoke.
NONLINEAR steel (B-H below). Excitation NI = 1000 / 3000 / 4500 / 5000 AT, to
sweep the saturation effect. Measurement points (mm): gap mid P1 (0,0,25.75),
gap edge P2 (12.5,5,25.75); average-Bz lines alpha-beta (center pole) and
gamma-delta (yoke). Quantity of interest = F_z on the center pole.

(An earlier note here said "two mu_r=1000 blocks 50x40x20mm" -- that was a
placeholder, NOT the real problem. Real TEAM-20 = pole+yoke+coil electromagnet
with NONLINEAR steel.) Exact Fig.1 dimensions are a low-res 1992 scan (labels
~50, 3.7, 50, 100, 98.5, 75, 39, gap ~1.x mm); take them from the lab TEAM-20
mesh/input before building an NGSolve model -- do NOT guess the flagship geometry.
Source PDF (owner-protected; decrypt with empty password via pymupdf
doc.save(..., encryption=PDF_ENCRYPT_NONE)):
W:\...\05_TEAM_benchmark\16_problem20\problem20_tech_3ax.pdf

## What it tests

- ★ Maxwell stress tensor implementation
- Force on iron in static B field
- Comparison of force methods (Maxwell, Coulomb, Kameari, VWP)
- Lab-implemented in `calc_em_force.py`

## Published reference (measured; 21 solutions from 11 groups)

  | NI [AT] | F_z [N] | Bz@P1 | Bz@P2 | <Bz> pole(a-b) | <Bz> yoke(g-d) |
  |---------|---------|-------|-------|----------------|----------------|
  | 1000    |   8.1   | 0.36  | 0.24  | 0.72           | 0.13           |
  | 3000    |  54.4   | 0.84  | 0.63  | 1.75           | 0.36           |
  | 4500    |  75.0   | 0.99  | 0.72  | 2.01           | 0.43           |
  | 5000    |  80.1   | 1.03  | 0.74  | 2.05           | 0.46           |

Steel B-H (Table 1, B[T]:H[A/m]) -- 0:0, 0.05:100, 0.10:153, 0.20:205, 0.30:233,
0.40:255, 0.50:285, 0.60:320, 0.70:355, 0.80:405, 0.90:470, 1.0:555, 1.1:673,
1.2:836, 1.3:1065, 1.35:1220, 1.4:1420, 1.45:1720, 1.5:2130, 1.55:2670, 1.6:3480,
1.65:4500, 1.7:5950, 1.75:7650, 1.8:10100, 1.85:13000, 1.9:15900, 1.95:21100,
2.0:26300, 2.05:32900, 2.1:42700, 2.15:61700, 2.2:84300, 2.25:110000, 2.3:135000.
The pole saturates to ~2 T (avg 2.05 T at 5000 AT) -- the >1.5 T tail dominates,
and the B-H must be extrapolated above 2.3 T (upper pole edge is very high).

## NGSolve / Radia recipe (validated lab approach -- REUSE this)

Geometry (netgen.occ, origin = yoke bottom centre, SI metres):
  yoke = Box(0,0,0)-(0.127,0.025,0.150) MINUS cavity X[0.025,0.102] Z[0.025,0.125]
         and pole-hole X[0.050,0.077] Z[0.125,0.150]   (both Y-through)
  pole = Box X[0.051,0.076] Y[0.0075,0.0175] Z[0.0265,0.125]  (1.5 mm bottom gap)
  coil = rectangular loop, inner 39x39 / outer 75x75 mm, wall 18 mm, Z[0.0267,0.1233],
         fillet R23 (outer) / R5 (inner) via MakeFillet on z-parallel edges.
Nonlinear steel: BSpline(2) over the 38-pt B-H -> energy density phi(B)=bh.Integrate();
  add it with Variation()/SymbolicEnergy -- a RAW add makes AssembleLinearization
  IGNORE the term (mu_r collapses, force wrong). Newton + line-search + load
  continuation (~12 NI steps). Coil J = per-wall (jx,jy) CF, J=NI/A_coil
  (A_coil = 18 x 96.6 mm). Open boundary: Kelvin transform (inner+outer sphere,
  Periodic Identify) OR a far Dirichlet box. Force = Maxwell stress on the
  'pole_surface' boundary with B on the AIR side (project curl A to L2(air), then
  BoundaryFromVolumeCF); Fz = -oint T_z ds(pole_surface).

## Cross-validation status (lab, Yano 2026-03)
Radia reproduces the measured force and SATURATES correctly:
  | NI [AT] | meas | Radia MST     | Radia Kelvin  |
  | 1000    | 8.1  | 6.86 (-15 %)  | 6.92 (-15 %)  |
  | 3000    | 54.4 | 47.5 (-13 %)  | 51.2 ( -6 %)  |
  | 4500    | 75.0 | 65.6 (-13 %)  | 71.0 ( -5 %)  |
  | 5000    | 80.1 | 70.1 (-12 %)  | 76.0 ( -5 %)  |
(Kelvin force method best, ~5 % at high NI.) The NGSolve A-Phi FEM is in progress
(NI=1000 -> 7.6 N vs 8.1). NOTE: one early A-method MST run over-shot to 210 N at
5000 AT -- that was a buggy / under-resolved run, NOT the converged answer; the
saturating steel must LEVEL the force off (measured 80.1 N). radia-ngsolve TEAM-20
flagship = finish the NGSolve A-form force at all NI (shipped
solve_magnetostatic_nonlinear + Maxwell-stress on the pole) and cross-validate
NGSolve vs Radia (-5..15 %) vs measured.

## radia-mcp pipeline result (2026-06-04) -- NI=1000 VALIDATED to 1.3 %
A radia-mcp port (Yano geometry, far-Dirichlet box, shipped solve_magnetostatic_Aform
+ air-side Maxwell stress on 'pole_surface') gives Fz = 8.20 N at NI=1000 (linear
mu_r=1000 in steel) vs measured 8.1 N = 1.3 %, BEATING Yano's own NGSolve A-Phi (7.6)
and Radia (6.86) at this unsaturated point. So the geometry + rectangular coil-J +
solve + force pipeline is validated. The pole body B came out 0.586 T (< 0.72
measured, because linear mu_r=1000 < the real ~1570), yet the FORCE matches: the net
pole force is a small difference of large opposing top/bottom/side surface stresses
dominated by gap fringing -- robust to body-flux error but sensitive to the surface
trace (hence the air-side L2 + BoundaryFromVolumeCF). The force normal/sign follow
Yano's -oint T_z ds('pole_surface'). NEXT: saturating cases (3000/4500/5000 AT) need
the nonlinear B-H + a FAST solver -- direct sparse-Cholesky is 185 s/solve (46k
order-2 tets) so Picard x ~30 iters ~ hours; use iterative CG+BDDC + Newton
to make nonlinear TEAM-20 tractable.

## radia-mcp NONLINEAR cross-validation (2026-06-04) -- COMPLETE, all NI <4.1 %
The shipped solve_magnetostatic_newton (Newton energy-min + CG/BDDC + load
continuation) on Yano's geometry + the 38-pt B-H, force = air-side Maxwell stress
on 'pole_surface', reproduces the MEASURED force at every NI (44k order-2 mesh):
  | NI [AT] | radia-mcp Fz | measured | err    | <|B|>_pole |
  | 1000    |  8.23 N      | 8.1      | +1.7 % | 0.59 T     |
  | 3000    | 56.62 N      | 54.4     | +4.1 % | 1.53 T     |
  | 4500    | 76.53 N      | 75.0     | +2.0 % | 1.77 T     |
  | 5000    | 81.49 N      | 80.1     | +1.7 % | 1.83 T     |
=> within +1.7..+4.1 % of measurement at ALL excitations -- BETTER than the lab's
own Radia (-5..-15 %) and NGSolve A-phi. The steel saturates correctly (pole
0.59 -> 1.83 T) so the force levels off (8 -> 57 -> 77 -> 81 N, not the ~200 N a
linear model gives). Newton converges in ~3-5 iters/load-step; CG/BDDC makes each
solve ~90-130 s (vs 185 s direct, and Picard would be hours). Shipped as the
validation case ``validate_team20_static_force`` (coarse NI=1000, ~15 s). This is
the FLAGSHIP: the radia-ngsolve magnetostatic force path is now validated end-to-end
against a measured, nonlinear, 3-D benchmark.

## Lab usage

★ THE canonical force benchmark for `radia_mcp.differential_forms.forces`.
All 7 implemented methods (Maxwell, Coulomb nodal, Kameari, Lorentz,
Arkkio, Henrotte, VWP) are validated here.

| Method | Lab implementation | Validates |
|--------|---------------------|-----------|
| Maxwell stress | `calc_em_force.py --method maxwell` | global force |
| Coulomb nodal | `calc_em_force.py --method coulomb_nodal` | per-node |
| Kameari (1993) | `calc_em_force.py --method kameari` | local distribution |
| Henrotte (2004) | `calc_em_force.py --method henrotte` | EM force density |

## Reference

[LOCAL] `05_TEAM_benchmark/problem20/problem20_3-D Static Force Problem.pdf`
[LOCAL] `05_TEAM_benchmark/problem20/problem20_tech_3ax.pdf` (technical spec)

## Cross-reference

- `radia_mcp.differential_forms.forces` — all 7 force methods
- `radia_mcp.electromagnet` — Hantila solver applied to problem 20
- `memory/feedback_em_force_team20_validation.md` (lab memory)
"""


PROBLEM_6 = r"""
# TEAM Problem 6 — Conducting Sphere in Uniform Sinusoidal Magnetic Field (★ analytical ref)

## Problem overview

A hollow conducting sphere (copper-like material, non-magnetic) placed in a uniform
sinusoidal magnetic field.  The induced eddy currents partially shield the interior
cavity.  Because the geometry has full spherical symmetry, an exact analytical solution
exists (modified spherical Bessel BVP), making this ideal for FEM validation.

```
Inner radius  R1 = 50 mm
Outer radius  R2 = 55 mm   (shell thickness 5 mm)
Conductivity  σ  = 5×10⁷ S/m  (copper-class)
Permeability  μᵣ = 1  (non-magnetic)
Applied field B₀ = 1 T  (uniform, z-directed)
Frequency     f  = 50 Hz  →  ω = 2π×50 rad/s
Skin depth    δ  = sqrt(2/(ωσμ₀)) ≈ 10.1 mm   (shell ≈ δ/2)
```

## Physics: scattered-field A-Φ formulation

Let A₀ = (B₀/2)(−y, x, 0) so that curl A₀ = B₀ẑ (background potential).
Solve for scattered Ã with Ã = 0 on an outer Dirichlet box:

```
∇×ν∇×Ã + jωσ(Ã + ∇Φ̃) = −jωσ A₀    (in shell only; σ=0 elsewhere)
jωσ(Ã + ∇Φ̃)·∇ψ          = −jωσ A₀·∇ψ  (gauge equation)
Ã = 0                                   (on outer box faces)
```

NGSolve `add_source` callback:
```python
A0 = CoefficientFunction((-B0*y/2, B0*x/2, 0.0))
def add_source(lf, test_fn):
    v_A, psi = test_fn
    lf += -1j * omega * sigma_cf * A0 * v_A       * dx
    lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx
```

Total B = curl(Ã) + B₀ẑ.  The ∇×ν∇×A₀ term vanishes on Dirichlet boundaries
(integration by parts), so no additional source from the background curl.

## Analytical reference (modified spherical Bessel BVP)

A_φ(r,θ) = f(r) sin θ, with B_z = 2f(r)/r (gives B_z = 2A uniform inside).

Regions:
- r < R1  (cavity):  f₁ = A·r
- R1<r<R2 (shell):   f₂ = C·i₁(κr) + D·k₁(κr),  κ = sqrt(jωσμ₀)
- r > R2  (outer):   f₃ = B₀r/2 + E/r²

where i₁(z) = cosh(z)/z − sinh(z)/z², k₁(z) = e^{−z}(1/z + 1/z²)

4×4 boundary-matching system (f, f' continuous at R1 and R2):
```python
kappa = np.sqrt(1j * omega * sigma * mu0)
z1, z2 = kappa*R1, kappa*R2

def i1(z):  return np.cosh(z)/z - np.sinh(z)/z**2
def k1(z):  return np.exp(-z)*(1/z + 1/z**2)
def i1p(z): return np.sinh(z)/z - 2*np.cosh(z)/z**2 + 2*np.sinh(z)/z**3
def k1p(z): return np.exp(-z)*(-1/z - 2/z**2 - 2/z**3)

M = [[R1, -i1(z1),       -k1(z1),       0      ],
     [1,  -k*i1p(z1),    -k*k1p(z1),    0      ],
     [0,   i1(z2),        k1(z2),       -1/R2**2],
     [0,   k*i1p(z2),     k*k1p(z2),    2/R2**3 ]]
rhs = [0, 0, B0*R2/2, B0/2]
A_coeff = np.linalg.solve(M, rhs)[0]
B_z_ref = 2 * A_coeff    # complex; |B_z_ref| ≈ 0.50 T, phase ≈ −60°
```

Thin-shell approximation (valid d ≪ δ): T ≈ 1/(1 + jωμ₀σ d R/3)
  → |T| ≈ 0.50, ∠T ≈ −60° for our parameters.

## NGSolve recipe (Dirichlet box, periodic=False)

```python
from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt

L = 0.25  # Dirichlet box half-length (truncation error < 1%)
sph_o = Sphere(Pnt(0,0,0), R2)
sph_i = Sphere(Pnt(0,0,0), R1)
box   = OrthoBrick(Pnt(-L,-L,-L), Pnt(L,L,L)).bc("outer")

shell = (sph_o - sph_i).mat("shell").maxh(0.003)  # 3 mm, ≥2 layers in δ/2
cav   = sph_i.mat("air_in")
air   = (box - sph_o).mat("air_out")

geo = CSGeometry()
geo.Add(shell); geo.Add(cav); geo.Add(air)
mesh = Mesh(geo.GenerateMesh(maxh=0.05, grading=0.4))
mesh.Curve(2)

nu       = CoefficientFunction([1/(4*pi*1e-7)] * len(mesh.GetMaterials()))
sigma_cf = mesh.MaterialCF({"shell": 5e7}, default=0.0)
A0       = CoefficientFunction((-B0*y/2, B0*x/2, 0.0))

gfAPhi = solve_eddy_current_harmonic_APhi(
    mesh, nu, sigma_cf, omega, add_source,
    order=2, precond="local", dirichlet="outer", periodic=False,
)
gfA = gfAPhi.components[0]
B_z_fem = complex(curl(gfA)[2](mesh(0.01, 0.0, 0.0))) + B0
```

## Regression test

`validation/force/validate_force_xval.py::validate_team6_sphere_eddy`
Pass criterion: |B_z_fem − B_z_ref| < 0.06 T (6 % of B₀)

## Reference values

| Quantity | Analytical (exact BVP, B₀=1T, 50Hz) |
|----------|---------------------------------------|
| B_z inside (complex) | 0.2329 − 0.4628j T |
| |B_z inside|          | 0.518 T |
| Phase shift           | −63.3° |

(Thin-shell approximation: |T| ≈ 0.50, phase ≈ −60° — close to exact.)

## Lab source paths

No Yano notebook for Problem 6.  Analytical solution self-contained in the test.
"""


PROBLEM_CYLINDER_AXIAL = r"""
# Conducting Cylinder in Uniform Axial Sinusoidal B Field (★ analytical ref)

## Overview

Analogous to TEAM Problem 2 (Infinitely Long Cylinder in a Sinusoidal Field).
A conducting non-magnetic cylinder (axis = z) immersed in a uniform axial B₀ẑ field.
Solved with the same scattered-field A-Φ approach as TEAM-6.

```
Cylinder: radius R = 50 mm, height H = 300 mm (finite, but 6R so mid-plane ≈ 2-D)
Material: σ = 1.6×10⁷ S/m (copper-like), μr = 1
Applied:  B₀ = 1 T, f = 50 Hz
Skin depth δ = √(2/(ωσμ₀)) ≈ 17.8 mm  →  R/δ ≈ 2.81 (partial penetration)
```

## Analytical solution (infinite cylinder)

B_z(r) = B₀ · I₀(κr) / I₀(κR)    where κ = √(jωσμ₀) = (1+j)/δ

- κR ≈ 2.81(1+j),  I₀(κR) ≈ -2.48 + 2.30j  → |I₀(κR)| ≈ 3.38
- B_z(0) = B₀/I₀(κR) ≈ -0.217 - 0.201j T,  |B_z(0)| ≈ 0.296 T
- ~70 % amplitude attenuation, ~-137 deg phase shift at centre  (R/delta = 2.81 → large delay)

I₀(z) = Σ_{k=0}^∞ (z/2)^{2k}/(k!)² — converges fast for |z| < 10.

## Scattered-field NGSolve recipe

Same source term as TEAM-6 sphere:
```python
A0 = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

def add_source(lf, test_fn):
    v_A, psi = test_fn
    lf += -1j * omega * sigma_cf * A0 * v_A    * dx
    lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx

gfAPhi = solve_eddy_current_harmonic_APhi(
    mesh, nu, sigma_cf, omega, add_source,
    order=2, precond="local", dirichlet="outer", periodic=False,
    tol=1e-8, maxsteps=2000, restart=200, reg=1e-10)

gfA     = gfAPhi.components[0]
B_z_fem = complex(curl(gfA)[2](mesh(0.001, 0.0, 0.0))) + B0
```

## Regression test

``test_cylinder_axial_eddy``  (order=2, cond.maxh=0.010 m, global maxh=0.05 m):
  Pass criterion: |B_z_fem(centre) − B_z_ref(centre)| < 0.05 T  (5 % of B₀)
"""


PROBLEM_SPHERE_SOLID_EDDY = r"""
# Solid Conducting Sphere in Uniform Axial Sinusoidal B Field (★ analytical ref)

## Overview

A solid (not hollow) conducting sphere in a uniform axial AC field — simpler 2-region BVP
than TEAM-6 (which uses a thin shell).  Analytical solution uses modified spherical Bessel
functions i₁(z) = cosh(z)/z − sinh(z)/z² and their real-number counterparts sinh/cosh.

```
Sphere:   radius R = 50 mm, solid conductor
Material: σ = 1.6×10⁷ S/m, μr = 1
Applied:  B₀ = 1 T, f = 50 Hz
Skin depth δ ≈ 17.8 mm  →  R/δ ≈ 2.81  (same parameters as cylinder test)
```

## Analytical solution

Inside the sphere A_total = α·i₁(κr)·sinθ eφ, regular modified spherical Bessel:
  i₁(z) = cosh(z)/z − sinh(z)/z²
Matching at r=R (continuity of Aφ and Bθ) gives:
  α = 3 B₀ R / (2 sinh(p))   where p = κR = (1+j) R/δ

B_z at sphere centre (r=0):
  As r→0: i₁(κr)→κr/3, so B_z(0) = 2ακ/3 = B₀ · p / sinh(p)

For p = 2.81(1+j):
  B_z_ref ≈ −0.209 − 0.432j T,  |B_z_ref| ≈ 0.480 T,  phase ≈ −116°

Limit checks: p→0 → B_z → B₀ (no shielding); p→∞ → B_z → 0 (perfect shielding). ✓

Note: unlike the cylinder, the sphere field is NOT uniform inside — only centre value given.

## Python formula

```python
kappa   = np.sqrt(1j * omega * sigma * mu0)
p       = kappa * R
B_z_ref = complex(B0 * p / np.sinh(p))
```

## NGSolve recipe

Same source term as TEAM-6 and cylinder:
```python
sigma_cf = mesh.MaterialCF({"conductor": sigma}, default=0.0)
A0       = CoefficientFunction((-B0 * y / 2, B0 * x / 2, 0.0))

def add_source(lf, test_fn):
    v_A, psi = test_fn
    lf += -1j * omega * sigma_cf * A0 * v_A       * dx
    lf += -1j * omega * sigma_cf * A0 * grad(psi) * dx
```

Geometry (CSG):
```python
sph  = Sphere(Pnt(0, 0, 0), R).mat("conductor").maxh(0.010)
air  = (OrthoBrick(Pnt(-L,-L,-L), Pnt(L,L,L)).bc("outer") - sph).mat("air")
```

Post-processing: B_z_fem = complex(curl(gfA)[2](mesh(0.001, 0.0, 0.0))) + B0

## Regression test

``test_conducting_sphere_axial_eddy``  (order=2, cond.maxh=0.006 m, global maxh=0.05 m):
  Pass criterion: |B_z_fem(centre) − B_z_ref(centre)| < 0.05 T  (5 % of B₀)
"""


PROBLEM_SPHERE_MAGNETOSTATIC = r"""
# Permeable Sphere in Uniform Magnetic Field — Scattered-Field Magnetostatics (★ analytical ref)

## Overview

A linearly permeable sphere (μr > 1, σ=0) in a uniform applied B₀ẑ field.
The exact (demagnetization factor N=1/3 for sphere) analytical result gives
a uniform field inside — ideal for validating the scattered-field A-form magnetostatic solver.

```
Sphere:   R = 50 mm
Material: μr = 10, σ = 0 (no eddy currents)
Applied:  B₀ = 1 T (z-directed, static)
```

## Analytical reference

B_z_inside = 3μr/(μr+2) × B₀  (exact, uniform inside sphere)
  For μr=10: B_z = 30/12 × 1 T = 2.5 T

This follows from N=1/3 demagnetization factor for a sphere:
  H_inside = H₀ − N × M_inside  →  uniform H inside  →  uniform B inside.

## NGSolve recipe (custom BilinearForm, NOT solve_eddy_current_harmonic_APhi)

```python
nu_cf   = mesh.MaterialCF({"sphere": nu0/mu_r}, default=nu0)  # ν = 1/μ
nu_pert = mesh.MaterialCF({"sphere": nu0*(1.0 - 1.0/mu_r)}, default=0.0)
B0_vec  = CoefficientFunction((0.0, 0.0, B0))

fes = HCurl(mesh, order=2, dirichlet="outer")
u, v = fes.TnT()
a = BilinearForm(fes)
a += nu_cf * InnerProduct(curl(u), curl(v)) * dx + 1e-8 * nu0 * InnerProduct(u, v) * dx
f = LinearForm(fes)
f += InnerProduct(nu_pert * B0_vec, curl(v)) * dx   # source: permeability jump
a.Assemble(); f.Assemble()
gfA = GridFunction(fes)
gfA.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * f.vec
```

Source derivation: scattered-field bilinear form.  Applied field B₀ = curl A₀ satisfies
∇×(ν₀∇×A₀) = 0 everywhere.  Inside sphere ν=ν₀/μr ≠ ν₀.  Subtract:
  ∇×(ν∇×Ã) = −∇×((ν−ν₀)∇×A₀) = +∇×(nu_pert × B₀)
The sign: nu_pert = ν₀(1 − 1/μr) > 0 for μr>1 — field is amplified.

Post-processing: B_z_fem = float(curl(gfA)[2](mesh(0.005, 0.0, 0.0))) + B0

## Regression test

``test_sphere_magnetostatic_scattered``  (order=2, sphere.maxh=0.008 m, global maxh=0.04 m):
  Pass criterion: |B_z_fem(centre) − 2.5 T| < 0.15 T  (6 % of B_z_inside)

Note: field inside a permeable sphere is *nearly* but not perfectly uniform due to FEM
discretisation.  Tolerance accounts for mesh-discretisation near sphere boundary.
"""


PROBLEM_21 = r"""
# TEAM Problem 21 P21a-0 — Tank-Wall Eddy Current (★ lab A-Φ Dirichlet-box test)

## Problem overview (TEAM P21 Family V.2009)

Models leakage-flux eddy currents in a large transformer tank wall.
Source: two rectangular stranded coils (HV/LV winding mockup), stacked in z.
Load:  non-magnetic steel plate to the side of the coil (+x direction).

```
Coil: outer=248 mm (R10 fillet), inner=200 mm (R45 fillet), height=217 mm each
      gap between coils = 24 mm,  turns N=300/coil,  f=50 Hz,  I=10 A rms
Case I: coil1 = +J (CCW, +z flux), coil2 = −J (CW, −z flux)
        → leakage flux loops through plate in x-direction (normal to plate face)
Case II: coil1 = coil2 = +J  → flux parallel to plate face

Plate P21a-0: 458 × 820 × 10 mm, 20Mn23Al non-magnetic steel
              σ=1.3889×10⁶ S/m,  μr=1
              skin depth δ = 1/√(π·f·μ₀·σ) ≈ 60 mm >> 10 mm plate → quasi-static
              (virtually no skin effect; coarse mesh adequate through thickness)
Plate position: 12 mm gap from coil face → x = 0.136 to 0.146 m
```

## Variants in the P21 family (V.2009 — 16 models)

| Model   | Steel type   | σ [S/m]   | μr  | Slits | δ [mm] |
|---------|-------------|-----------|-----|-------|--------|
| P21-B   | A3 magnetic | 6.484e6   | ~100| 0     | ~2.8   |
| P21a-0  | 20Mn23Al    | 1.3889e6  | 1   | 0     | ~60    |
| P21a-1  | 20Mn23Al    | 1.3889e6  | 1   | 1     | ~60    |
| P21b-MN | mixed       | —         | —   | 0     | —      |
| P21c-EM1| Cu shield   | 5.714e7   | 1   | —     | ~9.4   |

## NGSolve recipe — Dirichlet box (NO Kelvin transform)

Use `solve_eddy_current_harmonic_APhi` with `periodic=False`:

```python
# Material CFs  (all μr=1 for P21a-0 → ν=ν₀ everywhere)
nu       = CoefficientFunction([nu0] * len(mesh.GetMaterials()))
sigma_cf = mesh.MaterialCF({"plate": 1.3889e6}, default=0.0)

# Coil winding tangent: CCW around inner 200-mm square at origin
projx    = IfPos(x + 0.1, IfPos(x - 0.1,  0.1, x), -0.1)
projy    = IfPos(y + 0.1, IfPos(y - 0.1,  0.1, y), -0.1)
tau_coil = CF((projy - y, x - projx, 0)) / sqrt((projy-y)**2+(x-projx)**2+1e-30)

# Peak NI for 10 A rms: NI_peak = N * I_rms * √2
J_scale  = 300 * 10 * sqrt(2) / (0.024 * 0.217)   # coil_w=24mm, coil_h=217mm

def add_source(lf, test_fn):
    v, _ = test_fn
    lf += J_scale * tau_coil * v * dx("coil1")   # CCW → +z flux
    lf += (-J_scale) * tau_coil * v * dx("coil2") # CW  → −z flux  (Case I)

gfAPhi = solve_eddy_current_harmonic_APhi(
    mesh, nu, sigma_cf, omega, add_source,
    order=2, precond="local",
    dirichlet="outer",   # Dirichlet BC on outer air-box faces
    periodic=False,      # NO Kelvin/Periodic wrapping
    tol=1e-8, maxsteps=2000, restart=200, reg=1e-10
)
```

## Post-processing: eddy current loss

```python
gfA, gfPhi = gfAPhi.components
E_cf   = -1j * omega * (gfA + grad(gfPhi))   # E = −jω(A + ∇Φ) = J/σ
# Time-averaged power: P = ½ σ |E|²  integrated over conductor
P_eddy = (0.5 * Integrate(
    sigma_cf * InnerProduct(E_cf, Conj(E_cf)).real, mesh)).real
```

## Reference values (P21a-0, Case I, 10 A rms, 50 Hz)

| Quantity | Measured | Calculated | Note |
|---------|---------|-----------|------|
| Eddy loss P [W] | **9.17** | 9.31 | pure eddy (no hysteresis, μr=1) |
| Bx at plate face [×10⁻⁴ T] | — | — | x = 5.76 mm, y=0, vs z |

Other P21 family losses at 10 A rms Case I:
  P21-B (A3 magnetic): 11.97 W measured,  12.04 W calculated
  P21a-1 (1 slit):  3.40 W,  P21a-2 (2 slits): 1.68 W,  P21a-3 (4 slits): 1.25 W

## Regression test

``test_team21_eddy_loss`` (order=2, mesh_maxh=0.07 m, plate.faces.maxh=0.025 m):
  Pass criterion: |P_eddy − 9.17 W| < 3.5 W  (±38 % to tolerate coarse mesh)

## Lab source

Reference notebook (Yano, P21e variant — different coil/plate dimensions):
  W:\\00_CAE\\NGSolve\\矢野\\2026_03_30_TEAM_benchmark\\Problem21\\Problem_21_FEM_APhi.ipynb
TEAM P21 family spec (V.2009):
  W:\\03_文献・論文\\00_電磁界解析\\ベンチマーク問題(TEAM)\\problem21_Family (V.2009).pdf
"""


PROBLEM_21_P21B2N = r"""
# TEAM Problem 21 P21b-2N — Tank-Wall with Two Non-Magnetic Plates + Gap

## Geometry

Same coils as P21a family, but two narrower non-magnetic plates (20Mn23Al) separated by
a 200 mm gap centred on the coil axis:

```
Each plate: 458 × 248 × 10 mm  (same material, σ=1.3889×10⁶ S/m, μr=1)
Plate 1:  y from −348 mm to −100 mm  (outer edge to gap edge)
Plate 2:  y from +100 mm to +348 mm
Gap:      y from −100 mm to +100 mm  (200 mm wide, centred at y = 0)

Compare with P21a-0: single plate 458×820×10 mm (no gap, full width)
The gap reduces loss ~85 % relative to P21a-0 by allowing flux to pass through.
```

## NGSolve recipe (OCC geometry)

```python
half_gap, plate_w = 0.100, 0.248  # [m]

plate1 = Box(Pnt(plate_x0, -(half_gap+plate_w), 0.0), Pnt(plate_x1, -half_gap, plate_z))
plate2 = Box(Pnt(plate_x0,  half_gap,           0.0), Pnt(plate_x1, half_gap+plate_w, plate_z))
for p in (plate1, plate2):
    p.faces.maxh = 0.025
    p.mat("plate")

air_box = Box(
    Pnt(-half_o - margin, -(plate_y_max + margin), -margin),
    Pnt(plate_x1 + margin,  plate_y_max + margin,  plate_z + margin))
air_box.faces.name = "outer"
air = air_box - coil1 - coil2 - plate1 - plate2
air.mat("air")
```

## Reference values (P21b-2N, Case I, 10 A rms, 50 Hz)

| Quantity | Measured | Calculated | Note |
|---------|---------|-----------|------|
| Eddy loss P [W] | **1.38** | 1.37 | pure eddy (both plates non-magnetic, μr=1) |

## Regression test

``test_team21_p21b2n_eddy_loss`` (order=2, mesh_maxh=0.07 m, plate.faces.maxh=0.025 m):
  Pass criterion: |P_eddy − 1.38 W| < 0.6 W  (+-43 %)
"""


PROBLEM_21_P21A1 = r"""
# TEAM Problem 21 P21a-1 — Tank-Wall Eddy Current with 1 Central Slit

## Geometry

Identical coils and plate position to P21a-0, but the plate has a single
through-thickness horizontal slit cutting across most of its width:

```
Plate: 458 × 820 × 10 mm, 20Mn23Al (σ=1.3889×10⁶ S/m, μr=1)
Slit:  660 × 10 × 10 mm  (y_width × z_height × x_thickness)
       centred at y = 0,  z = 229 mm (mid-height of 458 mm plate)
       runs from y = −330 mm to y = +330 mm  (leaves 80 mm solid at each end)

Skin depth δ ≈ 60 mm >> 10 mm plate → quasi-static (thin-plate regime)
Slit blocks the dominant current circulation → ~63 % loss reduction vs P21a-0
```

## Physics: why the slit matters

In P21a-0 the main eddy-current loop circulates in the yz-plane through the
full plate height (458 mm).  The slit at mid-height (z ≈ 229 mm) cuts this path,
forcing current to detour through the narrow (10 mm) solid bridges at y = ±330 mm.
The resistance of those bridges is much higher → current drops → loss ~3.4 W vs 9.17 W.

## NGSolve recipe

Same solver and coil setup as P21a-0.  Only change: subtract a slit box from the plate.

```python
from netgen.occ import Box, Pnt, Glue, OCCGeometry

plate_x0, plate_x1 = 0.136, 0.146   # same as P21a-0
plate_y  = 0.410
plate_z  = 0.458
slit_y, slit_z = 0.330, 0.005        # half-extents: 660 mm wide, 10 mm tall
ctr_z = plate_z / 2                  # 0.229 m

plate_box = Box(Pnt(plate_x0, -plate_y, 0.0), Pnt(plate_x1, plate_y, plate_z))
slit_box  = Box(Pnt(plate_x0, -slit_y, ctr_z - slit_z),
                Pnt(plate_x1,  slit_y, ctr_z + slit_z))
plate = plate_box - slit_box
plate.faces.maxh = 0.015   # resolve slit edges (~10 mm height)
plate.mat("plate")
```

## Reference values (P21a-1, Case I, 10 A rms, 50 Hz)

| Quantity | Measured | Calculated | Note |
|---------|---------|-----------|------|
| Eddy loss P [W] | **3.40** | 3.34 | pure eddy (no hysteresis, μr=1) |

Other variants at same conditions:
  P21a-0 (no slit): 9.17 W,  P21a-2 (2 slits): 1.68 W,  P21a-3 (4 slits): 1.25 W

## Regression test

``test_team21_p21a1_eddy_loss`` (order=2, mesh_maxh=0.07 m, plate.faces.maxh=0.015 m):
  Pass criterion: |P_eddy − 3.40 W| < 1.5 W  (+-44 %)
"""


PROBLEM_13 = r"""
# TEAM Problem 13 — 3-D Nonlinear Magnetostatic Problem (★ lab core)

## Geometry (official — Nakata/Takahashi, COMPUMAG 1991)

A vertical iron sheet (3.2 x 50 x 126.4 mm) flanked by two C-shaped return-flux
sheets, driven by a rectangular coil (NI = 1000 or 3000 AT).  Half-model (z >= 0,
Neumann BC at z = 0: H_x = H_y = 0 at the symmetry plane).

```
Vertical sheet  : 3.2 x 50 x 126.4 mm, centred at origin (x = +/-1.6 mm)
C-sheets (top)  : two horizontal + two return arms around the coil
Coil            : 4 straight bricks (50 mm long) + 4 arc corners (r = 50 mm)
                  centred at (x = +/-50 mm, y = +/-50 mm), width 50 x 50 mm
Air gap (CRITICAL): 0.5 mm between vertical sheet edge (x = 1.6 mm) and C-sheet
                    (x = 2.1 mm) -- the dominant reluctance in the circuit.
```

CRITICAL MESHING NOTE: the 0.5 mm gap must be resolved (< 0.8 mm elements).
Without gap seeds, global maxh = 50 mm gives ~13k coarse elements → |B| ~ 15x too
low. Yano's `getMeshingParameters(maxh_global=100)` places hundreds of
`mp.RestrictH(P, r, 0.8)` seeds near the gap corners, giving a 37k-element mesh
that resolves the gap correctly.

## B-H curve (Kruger-Lehmann steel, 77 points, B: 0..5 T, H: 0..2.26e6 A/m)

B_KL = [0.0, 2.5e-3, 5e-3, ..., 5.0] (77 values)
H_KL = [0.0, 16, 30, ..., 2.26e6] (77 values)
energy_density = BSpline(2, [0]+B_KL, H_KL).Integrate()

## Coil current source

J0 = NI / A_coil = 1000 / 2500e-6 = 400,000 A/m^2  (A_coil = 50x50 mm^2)

Corner J (tangential): Jx = (y - yc) / r,  Jy = -(x - xc) / r
  corner centres: (xc, yc) = (+/-50mm, +/-50mm)
Brick J: constant per-face (e.g. brick_back: Jx = J0, Jy = 0)

Split J into two scalar CFs:
  J = J0 * CF_x(per_mat) * (1,0,0) + J0 * CF_y(per_mat) * (0,1,0)

This lets NGSolve evaluate the correct component per element.

## NGSolve recipe (validated -- REUSE this)

  fes = HCurl(mesh, order=2, dirichlet="outer")   # no nograds for our formulation
  energy_density = BSpline(2, [0]+B_KL, H_KL).Integrate()
  gfu = solve_magnetostatic_newton(
      mesh, J, energy_density, "iron",
      reg=1e-6, tol=1e-10, max_iter=15       # bddc, 8 iters -> converged
  )
  B = curl(gfu)

Half-model Neumann (natural) BC at z = 0 is automatic -- do NOT add Dirichlet
there. The A.n = 0 Dirichlet BC applies only to "outer" (far-field box).

## Reference values (NI = 1000 AT, measurements 1-7)

Average |Bz| in vertical sheet (x in [0,1.6mm], y in [-25,25mm]):
  z =  1 mm : 1.330 T
  z = 10 mm : 1.329 T
  z = 20 mm : 1.286 T
  z = 30 mm : 1.225 T
  z = 40 mm : 1.129 T
  z = 50 mm : 0.985 T
  z = 60 mm : 0.655 T

Measurement formula (use numerical integration on a grid):
  avg|Bz| = trapz_y( trapz_x( |B_z(x,y,z)| ) ) / (50mm * 1.6mm)

## Cross-validation status (radia-mcp, 2026-06-04, bddc reg=1e-6)

Solver: solve_magnetostatic_newton, bddc preconditioner, 37k elements (Yano mesh).
Newton: 8 iterations, solve time ~121 s on 12-core workstation.

  | z [mm] | radia-mcp | measured | err    |
  |--------|-----------|----------|--------|
  |   1    |  1.308 T  | 1.330 T  | -1.7 % |
  |  30    |  1.231 T  | 1.225 T  | +0.5 % |
  |  60    |  0.657 T  | 0.655 T  | +0.3 % |

avg |B| in iron: 0.9118 T. All within 2 % of measurement.
This matches Yano's own SparseCholeski result (1.317 / 1.233 / 0.664 T, 109 s)
to better than 0.7 %; the formulation is validated.

Shipped as validation case ``validate_team13_bflux``
(10x20 grid, tol=0.10 T, ~2 min, skips if Yano's geometry module unavailable).

## Lab source

Geometry + getMeshingParameters:
  W:\\00_CAE\\NGSolve\\矢野\\2026_03_30_TEAM_benchmark\\Problem13\\geometry.py
Measurement functions:
  W:\\00_CAE\\NGSolve\\矢野\\2026_03_30_TEAM_benchmark\\Problem13\\measurement.py
Reference notebook:
  W:\\00_CAE\\NGSolve\\矢野\\2026_03_30_TEAM_benchmark\\Problem13\\solution_of_team13_problem.ipynb
"""


PROBLEM_23 = r"""
# TEAM Problem 23 — Forces in Permanent Magnets (★ lab core)

## Geometry

Two permanent magnets with various pole arrangements; compute the
mutual force as a function of separation and orientation.

```
PM1: NdFeB, 30x30x10 mm, Br = 1.4 T
PM2: same, separation parametric
Pole arrangements: attractive, repulsive, parallel, perpendicular
```

## What it tests

- ★ PM force computation (no current source, just M)
- Treatment of constant-M material in force integral
- Verification of force on M-fixed body

## Lab usage

★ Canonical benchmark for PM force in `calc_em_force.py`.
Used to validate:
- `permanent_magnet_force.py` lab implementation
- Radia's `ObjHexahedron(verts, [Mx,My,Mz])` PM source
- Per-element force density export

## Reference

[LOCAL] `05_TEAM_benchmark/problem23/problem23_Forces in Permanent Magnets.pdf`

## Cross-reference

- `radia_mcp.differential_forms.forces.permanent_magnet_force` ★
- `radia_mcp.magnetic_materials.permanent_magnet` — PM datasheets
"""


PROBLEM_33B = r"""
# TEAM Problem 33b — Electric Local Force (★ lab per-node validation)

## Geometry

C-shaped electromagnet pulling a movable plunger. Iron + coil + plunger
with measured forces at multiple positions.

```
C-yoke: nonlinear steel (BH provided)
Coil: known turns + current
Plunger: same steel, movable along x
Position scan: x = 0..50 mm
```

## What it tests

- ★ Local force formulations (force per node / per element)
- Comparison between MULTIPLE published numerical methods
- "Experimental Validation of Electric Local Force Formulations"
  (V3 version is the most cited)
- Verification of PER-NODE Kameari force decomposition

## Lab usage

★ Used for `radia_mcp.differential_forms.forces.local_methods`:
- Per-node Kameari (Kameari 1993)
- Iino-Okamoto two-stage error correction
- Eggshell method
- Sensitivity-method VWP (Pile 2018)

These methods provide DISTRIBUTION of force, not just global integral.

## Reference

[LOCAL] `05_TEAM_benchmark/problem33b/problem-33b_V3_Experimental Validation of Electric Local Force Formulations.pdf`

## Cross-reference

- `radia_mcp.differential_forms.forces.local_methods` ★
- `radia_mcp.differential_forms.forces.iino_okamoto`
- Lab memory: `tests/panels/test_calc_em_force_team33b_golden.py`
"""


PROBLEM_28 = r"""
# TEAM Problem 28 — Electrodynamic Levitation Device

## Geometry

Conducting plate above an axisymmetric coil set with AC excitation.
The plate floats due to eddy current Lorentz force.

```
Plate: aluminum, 50x50x5 mm
Coil: 3-axis axisymmetric (multiple windings)
Excitation: AC at 50 Hz, 10 A rms
Equilibrium height: ~5 mm
```

## What it tests

- ★ Levitation force balance (gravity vs eddy Lorentz)
- 3D harmonic eddy current with motion (force depends on height)
- Comparison with measured levitation height

## Lab usage

- Cross-link to `radia_mcp.maglev` (future MCP) for maglev physics
- Cross-link to `radia_mcp.ih` for eddy current in similar config

## Reference

[LOCAL] `05_TEAM_benchmark/problem28/` (11 files, 15 MB)
- problem28 description
- numerical solutions, measured data

## Cross-references

- `radia_mcp.differential_forms.forces.lorentz` — Lorentz force
- `radia_mcp.ih.sibc` — analogous eddy + force physics
"""


def get_force_motion_knowledge(topic: str = "overview") -> str:
    """Dispatch TEAM force / motion / nonlinear-B / eddy-current problem topics.

    Topics:
        overview             - All problems summary (DEFAULT)
        problem_7            - ★ 3-D Eddy Current A-Phi (Kelvin, lab A-Phi test)
        cylinder_axial       - ★ Conducting Cylinder in Axial B (TEAM-2 analog, I₀ Bessel ref)
        sphere_solid_eddy    - ★ Solid Conducting Sphere in Axial B (sinh/cosh Bessel ref)
        sphere_magnetostatic - ★ Permeable Sphere in Uniform B (demagnetization factor ref)
        problem_21           - ★ Tank-Wall Eddy Current P21a-0 (Dirichlet box)
        problem_21a1         - ★ P21a-1: same + 1 central slit (3.40 W ref)
        problem_21b2n        - ★ P21b-2N: two non-magnetic plates with gap (1.38 W ref)
        problem_13           - ★ 3-D Nonlinear C-Core B-Flux (lab core)
        problem_20           - ★ 3-D Static Force Problem (lab core)
        problem_23           - ★ Forces in Permanent Magnets (lab core)
        problem_33b          - ★ Electric Local Force (per-node validation)
        problem_28           - Electrodynamic Levitation
        all                  - Everything
    """
    topic = topic.lower().strip().replace("-", "_")
    if topic in ("overview", "all_force"):
        return OVERVIEW
    if topic in ("problem_6", "6", "problem6", "sphere_shell_eddy", "sphere_eddy", "scattered_field"):
        return PROBLEM_6
    if topic in ("cylinder_axial", "cylinder", "problem_2", "2", "problem2", "team2"):
        return PROBLEM_CYLINDER_AXIAL
    if topic in ("sphere_solid_eddy", "6s", "solid_sphere", "sphere_solid"):
        return PROBLEM_SPHERE_SOLID_EDDY
    if topic in ("sphere_magnetostatic", "6m", "permeable_sphere", "demagnetization"):
        return PROBLEM_SPHERE_MAGNETOSTATIC
    if topic in ("problem_7", "7", "problem7", "eddy", "eddy_current", "aphi"):
        return PROBLEM_7
    if topic in ("problem_21", "21", "problem21", "p21a", "tank_wall", "dirichlet_box"):
        return PROBLEM_21
    if topic in ("problem_21a1", "21a1", "p21a1", "p21a_1", "slit_1", "one_slit"):
        return PROBLEM_21_P21A1
    if topic in ("problem_21b2n", "21b2n", "p21b2n", "two_plates", "gapped_plate"):
        return PROBLEM_21_P21B2N
    if topic in ("problem_13", "13", "problem13", "ccore", "bflux"):
        return PROBLEM_13
    if topic in ("problem_20", "20", "problem20", "static_force"):
        return PROBLEM_20
    if topic in ("problem_23", "23", "problem23", "pm_force"):
        return PROBLEM_23
    if topic in ("problem_33b", "33b", "problem33b", "local_force"):
        return PROBLEM_33B
    if topic in ("problem_28", "28", "problem28", "levitation"):
        return PROBLEM_28
    if topic == "all":
        return "\n\n".join([OVERVIEW, PROBLEM_6, PROBLEM_CYLINDER_AXIAL,
                            PROBLEM_SPHERE_SOLID_EDDY, PROBLEM_SPHERE_MAGNETOSTATIC,
                            PROBLEM_7, PROBLEM_21, PROBLEM_21_P21A1, PROBLEM_21_P21B2N,
                            PROBLEM_13, PROBLEM_20, PROBLEM_23, PROBLEM_33B, PROBLEM_28])
    return (f"Unknown topic '{topic}'. Available: overview, problem_6, cylinder_axial, "
            "sphere_solid_eddy, sphere_magnetostatic, problem_7, "
            "problem_21, problem_21a1, problem_21b2n, problem_13, "
            "problem_20, problem_23, problem_33b, problem_28, all.")
