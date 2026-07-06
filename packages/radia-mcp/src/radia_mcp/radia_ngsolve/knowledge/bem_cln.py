"""BEM-CLN knowledge base for Radia MCP server.

Per-element multipole CLN with Schur-F termination, assembled into
multi-conductor systems via integral-equation Green's function coupling.

Backs the IEEE Transactions on Magnetics submission "Universal Cauer-SIBC
Composition via Schur Complement" (Sugahara, Nagamine, Hane, 2026)
sections V.G (DOF accounting) and V.H (verification).

The framework extends single-conductor Schur-F augmentation to multi-
conductor systems with:
  - Polarizability alpha(s) = volume - Y_cln(s) / sigma as the natural
    per-element variable (DC limit alpha=0, PEC limit alpha=V built-in)
  - 1/D^2 (2D) or mu_0/(4 pi D^3) (3D) inter-element Green's function
  - No phenomenological saturation factor (alpha is bounded by construction)

Companion files:
  - docs/cln/BEM_CLN.md  (full theoretical exposition)
  - validation_test/maglev/research_cln/bem_cln_*.wls  (Mathematica verification)
"""

BEM_CLN_OVERVIEW = """
# BEM-CLN: Per-Element Multipole CLN with Schur-F Termination

## What this is

A multi-element extension of the single-conductor Schur-F augmented
Cauer Ladder Network (CLN), suitable for assembling fast eddy-current
reduced-order models (ROMs) of small clusters of conductors using
natural open-boundary integral-equation coupling.

For a system of N conductors:

  Per-element DOF = N_Cauer + 1   (Cauer rungs + 1 Schur-F sqrt(s) block)
  Total DOF       = N (N_Cauer + 1)
  Coupling matrix = N x N (dense within element list, sparse overall)

Compare with FEM-element CLN (Takahashi 2024): ~60 N_Cauer DOFs per element
and bounded-domain only.  Compare with high-order curved FEM (p=3): ~60+
DOFs per element plus Kelvin/PML for open boundary.  BEM-CLN delivers
natural open boundary + exact SIBC asymptote at modest DOF cost.

## Mathematical structure

### Per-element building block (whole-conductor CLN + Schur-F)

For a single isolated conductor in applied transverse B_0:

  Y_cyl(s) = pi a^2 sigma * 2 I_1(gamma a) / (gamma a I_0(gamma a))     (2D cyl)
  Y_cuboid(s) = sigma * Sum_{m,n,p odd} g^2_mnp / (1 + s tau_mnp)        (3D cuboid)
  Y_cyl_S(s) = Y_CLN_N(s) + K_SIBC^(cyl) sqrt(s) / (s + d)              (Schur-F)

### Polarizability (the right BEM variable)

  alpha(s) = "conductor volume" - Y_cln(s) / sigma   (per-element)
           = a^2 - Y_cyl(s)/(pi sigma)   (2D)
           = V   - Y_cuboid(s)/sigma     (3D)

Limits built in:
  alpha(0)        = 0           (no eddy current at DC)
  alpha(infinity) = a^2 or V    (PEC diamagnetic shielding)
  |alpha(s)|     <= a^2 or V   (bounded -> no saturation factor needed)

### Assembly (N elements)

  c_n = alpha(s) * B^local_n       (induced dipole strength)
  B^local_n = B_0 + sum_{m != n} (Green's function kernel) * c_m

In matrix form:
  (I -/+ alpha(s) G) c = alpha(s) B_0 1_vec

with:
  2D coupling:  G_ij = 1 / D_ij^2          (sign: -)
  3D coupling:  G_ij = mu_0/(4 pi D_ij^3)  (sign: +, destructive 3D dipole)

System size: N x N.  Total polarizability: alpha_total = sum_n c_n / B_0.

## Limits

  DC:                alpha_total(0, D) = 0
  Far field:         alpha_total(s, D -> infty) = N alpha_single(s)
  PEC + coupling:    alpha_total(s -> infty, D) = N V / (1 +/- V * lambda_max(G))
  SIBC tail:         The K_SIBC sqrt(s) block in each per-element Schur-F
                     gives the 1/sqrt(s) approach to the PEC limit
                     structurally; the assembled system inherits this
                     asymptote automatically.

## Phase 1 -> 2 -> 2.5 -> 3 progression (May 2026)

  Phase 1   (2D 2-cyl, phenomenological saturating coupling) -> superseded
  Phase 2   (2D N-cyl scaling, same phenomenological) -> superseded
  Phase 2.5 (2D 2-cyl rigorous, polarizability) -> canonical for 2D
  Phase 3 A (3D quadrupole, Theorem 3 verification) -> Mathematica done
  Phase 3 B (3D 2-cuboid BEM-CLN with 3D dipole coupling) -> Mathematica done

  Outstanding: NGSolve A-formulation full-FEM cross-check (Phase 3 B+
  reference), real induction-heating workpiece N >= 5 multi-element
  example, real-time Verilog-A SPICE realisation of the sqrt(s) block.

## Key result

The Schur block construction is *structurally invariant* across
single-conductor (V.D vector decomposition), inhomogeneous applied field
(V.E Theorem 3), and multi-conductor BEM assembly (V.G, V.H): only the
scalar K_SIBC (or polarizability tensor) rescales.  This is what makes
the framework a candidate BEM-CLN solver kernel rather than just a
single-conductor ROM technique.
"""

BEM_CLN_2D_RIGOROUS = """
# Phase 2.5 Rigorous 2D BEM-CLN (canonical)

## Setup

Two parallel infinite Cu cylinders, radius a, separated by D,
applied uniform transverse B_0 in y direction.

  alpha(s) = a^2 [1 - 2 I_1(gamma a) / (gamma a I_0(gamma a))]
          = a^2 - Y_cyl(s)/(pi sigma)
  alpha(0) = 0,  alpha(infinity) = a^2

## 2D dipole-dipole coupling

A 2D dipole c cos(theta)/r at origin produces:

  Delta B_y(D, 0) = c / D^2     (perpendicular to applied field direction)

## 2-cylinder system (by symmetry m_1 = m_2 = m)

  m (1 - alpha(s)/D^2) = alpha(s) B_0
  alpha_total = 2 alpha(s) / (1 - alpha(s)/D^2)

  Stability: alpha/D^2 < 1 for D > a (bounded by alpha <= a^2)

## N-cylinder linear chain (D_nm = |n-m| D)

  (I - alpha(s) G_2D) c = alpha(s) B_0 1
  (G_2D)_ij = 1 / (|i-j| D)^2 for i != j

  Solve N x N system at each frequency s, compute alpha_total = sum c.

## Verified properties (5 x 2 mm Cu, basis-4 Schur-F per cylinder)

  - DC:           alpha_total(0) = 3 x 10^-16        (machine zero)
  - Far field:    alpha_total/(N alpha) = 1.0 + 10^-15  at D = 200 m
  - Stability:    a^2 lambda_max(G) = 0.235 < 1      (well-conditioned)
  - PEC + coupling: alpha_total -> N a^2 (1 + coupling enhancement)
  - Schur/Bessel ratio at f=10^6 Hz: 0.998-0.999  across N=2,3,5,10
    (consistent 0.13% undershoot, no N-amplification)

The rigorous 1/D^2 coupling on alpha (which has DC limit alpha=0 built
in) requires no phenomenological saturation factor: the earlier Phase 1
g(s,D) = mu_0/(2 pi D^2) * s/(s + s_sat) is superseded.
"""

BEM_CLN_3D = """
# Phase 3 B Rigorous 3D BEM-CLN (Mathematica framework)

## Setup

Two finite Cu cuboids of dimensions a x b x c, separated by D in some
direction, applied uniform transverse B_0.

  alpha_3D(s) = V - Y_cuboid(s) / sigma     (V = a b c)
  alpha_3D(0) = 0
  alpha_3D(infinity) = V (PEC volume exclusion)

## 3D dipole-dipole coupling

A point dipole m_y at origin produces, at distance D along x:

  Delta B_y(D, 0, 0) = -mu_0 m_y / (4 pi D^3)

The minus sign comes from the (3 cos^2(theta) - 1) angular factor of
the 3D dipole field for moments perpendicular to the separation
direction; this is *opposite* to the 2D case where coupling is
positive (constructive).

## N-cuboid system

  (I + alpha(s) G_3D) c = alpha(s) B_0 1
  (G_3D)_ij = mu_0 / (4 pi |r_i - r_j|^3)

  Note: + sign in front of alpha G (vs - in 2D), reflecting 3D dipole
  geometry.

## Practical magnitudes

For typical 3D BEM-CLN parameters:
  Cuboid 5 x 2 x 1 mm Cu (V = 1e-8 m^3)
  D = 15 mm separation
  V * mu_0 / (4 pi D^3) = 1e-8 * 1e-7 / 3.4e-6 ~ 3e-10

So 3D BEM-CLN coupling magnitude at typical industrial separations is
~10^-10 (negligible relative to alpha ~ V ~ 1e-8).  The framework is
*structurally* valid in 3D, but the practical coupling contribution
becomes visible only at:
  - Very close conductors (D approaching conductor size)
  - Larger conductors (cm-scale or above)
  - High-multipole-order applied fields (where coupling structure
    differs from simple dipole-dipole)

This is consistent with the physical observation that in 3D, parallel
conductors at typical spacing affect each other primarily through their
shared volume of curl-A (mutual inductance), not through dipole-dipole
external coupling.

## Verified

  - DC limit: alpha_total(s->0) on 2-cuboid = numerical zero
  - Far-field: alpha_total/(N alpha_single) = 1.0 at D = 100 m
  - PEC limit: alpha_total/(N alpha_single) = 0.9999... at D = 15 mm
    (coupling correction is ~3e-10, negligible at this scale)
  - Foster vs Schur-F per-cuboid building blocks: agreement to 0.99 in
    wall band, 0.9999 at f = 10^12 Hz

## Outstanding

  - NGSolve A-formulation full-FEM reference solve to cross-check the
    3D coupling derivation (Phase 3 B+ task)
  - Real industrial multi-element example (induction-heating workpiece
    + multi-turn coil with N >= 5 elements)
  - Inclusion of quadrupole-order per-element response (the leading
    correction beyond uniform applied B, per Theorem 3 of the manuscript)
"""

BEM_CLN_NOTEBOOK_INDEX = """
# BEM-CLN reference scripts (Mathematica)

Path: validation_test/maglev/research_cln/

  bem_cln_2cylinder.wls               (Phase 1, legacy: 2D phenomenological coupling)
  bem_cln_Ncylinder.wls               (Phase 2, legacy: N-cyl chain phenomenological)
  bem_cln_2cylinder_rigorous.wls      (Phase 2.5 canonical: 2D rigorous polarizability)
  bem_cln_2cuboid_3D.wls              (Phase 3 B: 3D 2-cuboid with mu_0/(4 pi D^3))
  cuboid3D_quadrupole.wls             (Theorem 3 2D: u_ext = (y^2 - x^2)/2 on long cuboid)
  cuboid3D_quadrupole_3D.wls          (Theorem 3 3D: u_ext = xy on 5x2x1 mm cuboid)
  cuboid3D_vector_decomposition.wls   (Theorem 2: polarisation-resolved K_SIBC)

CSV outputs in supplement/data/.

Reference manuscript:
  Sugahara, Nagamine, Hane (2026), IEEE Trans. Magn., submitted.
  Sections V.G (DOF accounting), V.H (BEM-CLN verification).
"""


def get_bem_cln_documentation() -> str:
    """Return the full BEM-CLN documentation."""
    return "\n\n".join([
        BEM_CLN_OVERVIEW,
        BEM_CLN_2D_RIGOROUS,
        BEM_CLN_3D,
        BEM_CLN_NOTEBOOK_INDEX,
    ])
