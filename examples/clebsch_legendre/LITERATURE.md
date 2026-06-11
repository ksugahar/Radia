# Prior-Art Survey: Inverse-Variable (Hodograph/Legendre) Formulations of Clebsch Magnetostatics

Literature assessment (web survey, 2026-06-11) for the two formulations in
this directory:

- **Variant 1 (full swap)**: coordinates `(phi, psi, z)`, unknowns
  `x(phi,psi,z)`, `y(phi,psi,z)`, energy
  `W = Int (1+x_z^2+y_z^2)/(2 mu0 J) dphi dpsi dz`.
- **Variant 2 (partial swap)**: coordinates `(x, y, psi)`, unknowns
  `z(x,y,psi)`, `phi(x,y,psi)`; 2D limit
  `d/dx(2 z_x/z_psi) = d/dpsi((1+z_x^2)/z_psi^2)`.

## Verdict

The **mechanism** (swap dependent/independent variables; flux/stream labels
as coordinates; geometry as unknown; energy minimization) is **classical in
several neighboring fields**. The **combination** here — raw Clebsch pair +
one Cartesian coordinate in a slab (non-toroidal) setting, the specific
vacuum energy functional, the partial `z <-> psi` swap with the in-surface
potential as second unknown, and the application to current-free
magnet-gap / pole-face design — was **not found in any of the five
literature directions searched**. Novelty should be claimed at the level of
the combination/application, never at the level of the mechanism.

## What is established (cite this lineage)

### 2D potential-plane methods (full swap) — classical and standard

| Work | What it does |
|---|---|
| Rogowski, Arch. Elektrotech. 12, 1 (1923) | Electrode/pole = chosen equipotential of a prescribed field |
| Weber, *Electromagnetic Fields* (1950); Binns & Lawrenson (1963/1973) | `z = f(w)`, `w = phi + i psi` plane as the rectangle where every 2D problem is trivial (curvilinear squares) |
| Beth, J. Appl. Phys. 37, 2568 (1966) | Complex/analytic representation of 2D magnet fields |
| Halbach, Nucl. Instr. Meth. 64, 278 (1968); 74 (1969); 78, 185 (1970) | Conformal-mapping evaluation/design of accelerator magnet poles; still in use (e.g. dipole shim design, NIM-A 2019, arXiv:1801.05470) |
| Chaplygin (1902); Courant & Friedrichs (1948); Bers (1958) | Hodograph transformation, textbook |
| Henrot & Pierre, Quart. Appl. Math. 49 (1991) | 2D free-boundary electromagnetic shaping via conformal map (geometry unknown) |

### Partial swap — classical in FLUID mechanics only

| Work | What it does |
|---|---|
| **von Mises, ZAMM 7, 425 (1927)** | Boundary layer in `(x, psi)` coordinates, `y(x,psi)` unknown — the exact 2D twin of Variant 2 |
| Martin, Arch. Rat. Mech. Anal. 41, 266 (1971) | Viscous flow in `(phi, psi)` curvilinear coordinates |
| Huang & Dulikravich, CMAME 59, 155 (1986); AIAA 91-0189 | "Stream-function-coordinate" (SFC) formulation; the SFC equation is structurally our 2D-limit equation |
| Butterweck & Pozorski, Acta Mech. 224, 1801 (2013) | Viscous inverse design of wall shape `y(x,psi)` |
| Stanitz, NACA TR-1115 (1953); NASA CR-3288 (1980) | Inverse channel design in `(phi, psi)` (2D) and potential + two stream functions (3D) — geometry recovered from prescribed wall velocity |
| Greywall, CMAME (1983), JCP (1988), ASME JFE 115, 233 (1993) | 3D flows with one streamwise spatial coordinate + two stream-surface labels as independents — the same independent-variable structure as Variant 1, for flows, no variational form |

No magnetostatic transplant of the partial swap was found.

### Plasma physics: inverse flux coordinates + energy minimization — classical

| Work | What it does |
|---|---|
| Zakharov & Shafranov, Sov. Phys. Tech. Phys. 18, 151 (1973); Rev. Plasma Phys. 11 (1986) | Axisymmetric inverse variables `r(a,theta)`, `z(a,theta)` |
| Bauer, Betancourt & Garabedian, *A Computational Method in Plasma Physics* (1978) | Variational inverse-coordinate equilibrium (BETA); includes vacuum-shell chapters |
| DeLucia, Jardin & Todd, JCP 37, 183 (1980); Lao, Hirshman & Wieland, Phys. Fluids 24, 1431 (1981) | Inverse-variable / variational moment Grad-Shafranov |
| **Hirshman & Whitson, Phys. Fluids 26, 3553 (1983)** (VMEC) | "Inverse coordinate representation `x = x(rho,theta,zeta)`", minimizes `W = Int (B^2/2mu0 + p/(gamma-1)) dV` by steepest descent. At `p = 0` this is our Variant-1 idea in toroidal flux coordinates. Modern: DESC (Dudt & Kolemen 2020) |
| D'haeseleer, Hitchon, Callen & Shohet, *Flux Coordinates and Magnetic Field Structure* (1991) | Clebsch `(alpha, beta, l)` as curvilinear coordinates — descriptive only (no position-as-unknown solve). The canonical reference for the coordinate system itself |
| Cary & Littlejohn, Ann. Phys. 151, 1 (1983) | Field-line flow is Hamiltonian with the Clebsch pair canonical and the third coordinate as time — Variant 1 is a transformation to these canonical variables |

### Vacuum fields with position-as-unknown (closest single line)

| Work | What it does |
|---|---|
| **Boozer, Phys. Plasmas 26, 102504 (2019)** "Curl-free magnetic fields for stellarator optimization" (arXiv:1906.06807) | Curl-free annular field constructed in Boozer coordinates with position `x(psi,theta,phi)` unknown; curl-free + div-free imposed by equating two representations of B. **Closest prior art.** Toroidal, Fourier, flux/angle variables (not the raw Clebsch pair), no slab/pole-design setting |
| Boozer, Phys. Plasmas 26, 042104 (2019) (arXiv:1812.05673) | Clebsch `(alpha, beta, l)` coordinates with Cartesian `x(alpha,beta,l,t)` dependent — kinematic (ideal evolution), not a magnetostatic BVP |
| Garren & Boozer, Phys. Fluids B 3, 2822 (1991); Landreman & Sengupta, JPP 84 (2018); + Plunk JPP 85 (2019) | Near-axis direct construction: prescribed `B`(Boozer coords) -> surface shapes; vacuum case standard |
| Giuliani, Wechsung, Stadler, Cerfon & Landreman, JPP 88 (2022) (arXiv:2203.03753) | Numerical computation of vacuum-field surfaces parametrized in Boozer coordinates (residual solve) |
| Sengupta & Weitzner, JPP (2019); Cary, PRL 49, 276 (1982); Grad (1967) | Existence theory: vacuum fields with nested toroidal flux surfaces are generically obstructed (Grad's conjecture). Caveat literature — the slab/boundary-driven magnet-gap setting differs |
| Cally, JCP 89, 388 (1991) | 2D free-boundary magnetohydrostatics in inverse flux coordinates (multigrid) — closest published 2D relative of our 2D-limit equation |

### Magnet/pole engineering literature

All found pole-design work is **optimization over parametrized geometry**
(adjoint shape sensitivity: Park, Coulomb & Hahn IEEE Trans. Magn. 27/29;
ROXIE/Russenschuck "inverse field computation" = harmonic fitting; MRI
pole/shim design = regularized linear inversion). **No true
inverse-coordinate PDE formulation** (geometry as the unknown of a
flux-coordinate PDE) exists in this literature. Japanese-language searches
(磁極形状 逆問題, Clebsch ポテンシャル 静磁場, 磁束座標) also returned
nothing of the kind.

## What appears novel (as of this survey)

1. **Variant 2 in 3D**: the partial `z <-> psi` swap with unknowns
   `z(x,y,psi)` + in-surface potential `phi(x,y,psi)` has **no published
   counterpart in any field** (its 2D restriction is von Mises / SFC, for
   flows).
2. **Application to pole-face / magnet-gap design** (pole = prescribed flux
   surface, geometry solved as PDE unknown): unpublished; the 2D
   conformal-mapping pole tradition (Rogowski/Halbach) is the historical
   special case being generalized.
3. **Variant 1's specific form**: raw Clebsch pair + Cartesian `z` (slab
   topology) with the vacuum energy functional
   `Int (1+x_z^2+y_z^2)/(2 mu0 J)`. The idea is the `p = 0` limit of
   VMEC/BETA transplanted out of the toroidal angle representation; that
   transplant is unclaimed.

## Verification caveats (do before citing in print)

- Boozer PoP 26, 102504 (2019), Giuliani et al. (2022), Sengupta et al.
  (2021): equation-level overlap was assessed from search extracts only
  (publisher/arXiv fetches were blocked, HTTP 403, in the survey
  environment). Pull the PDFs and confirm before any written novelty claim.
- Stanitz NASA CR-3288 (1980): exact 3D dependent-variable set reported
  from secondary sources; verify against the report PDF.
- Halbach NIM 78, 185 (1970) details confirmed only via secondary
  citations.
