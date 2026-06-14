r"""Why open-boundary methods stay accurate on COARSE meshes -- a spectral
(DtN-matrix) explanation of the Kelvin-transformation coarse-mesh accuracy.

Kameari's classic demonstration shows, by mesh REFINEMENT, that the Kelvin
transformation already gives good accuracy on relatively coarse meshes.  This
module reframes that *empirical* observation as a *property of the DtN matrix*
measured across mesh sizes:

  * Every open-boundary closure targets ONE continuous operator on the
    truncation surface Γ -- the exterior Dirichlet-to-Neumann (Steklov-Poincaré)
    map Λ_ext.  Kelvin transform and BEM coupling are DISCRETIZATIONS of Λ_ext
    (exact in the continuum, h→0); PML / infinite elements / asymptotic Robin BC
    are MODEL surrogates with a fixed (mesh-independent) modal error.
  * On a sphere of radius R, the CONTINUOUS Λ_ext is diagonalised by the
    spherical harmonics with the mesh-independent eigenvalue ladder
    λ_n = −(n+1)/R  (degree n, multiplicity 2n+1).  (This ladder is the
    yardstick; the DISCRETE matrix Λ_h does not share its mesh-independence.)
  * The discrete DtN matrix Λ_h reproduces the LOW-degree eigenvalues accurately
    already on the coarsest mesh (from a tiny floor that still falls steeply,
    ~O(h⁴), under refinement); the per-mode error grows with degree n.  Because
    the exterior field of a compact source is dominated by the low multipoles,
    the coarse mesh already resolves everything that matters -- which matches the
    accuracy Kameari saw, now stated as a spectral fact rather than a refinement
    curve.  Measured BOTH ways: the BEM Λ_h spectrum directly, AND the Kelvin
    closure's effective DtN (volume FEM on the inverted ball) -- two
    discretisations of the one Λ_ext, both coarse-mesh accurate for the low modes.

Companion code:
  radia_ngsolve.bem_integral (BEM side -- the boundary operator spectrum)
    - laplace_exterior_dtn():     assemble Λ_h = V⁻¹(−½M + K)
    - exterior_dtn_spectrum():    eigenvalues of Λ_h matched to −(n+1)/R
    - dtn_spectrum_vs_mesh():     per-degree eigenvalue error vs mesh size
  radia_ngsolve.fem_bem_coupling (Kelvin side -- the volume-FEM realisation)
    - kelvin_dtn_eigenvalue():    effective DtN of the Kelvin closure for mode n,
                                  order>=n kills the poly error -> geometry floor
                                  (~5-6 digits in 3D, deeper in 2D);
                                  dim=3 (sphere, −(n+1)/R) or dim=2 (circle, −n/R,
                                  the static-apparatus / rotating-machine cross-section)
    - kelvin_vs_exact_open_bc_error(): the Kelvin open-BC error ISOLATED from the
                                  interior FEM error (shared mesh, swap only Γ op)
    - kelvin_openbc_error_vs_exterior_mesh(): the isolated open-BC error vs the
                                  EXTERIOR (Kelvin ball) mesh size, interior fixed
    - kelvin_twosphere_shell_dipole(): the lab's REAL two-sphere periodic Kelvin
                                  (validation: reproduces the exact dipole)

Related knowledge: kelvin_transformation (Kelvin material modulation),
fem_bem_schur (using Λ_ext as an exact open BC).
"""

DTN_COARSE_MESH_OVERVIEW = r"""
# Coarse-Mesh Accuracy of Open BCs = a Property of the DtN Matrix

## The question (Kameari, reframed)

Kameari demonstrated the Kelvin transformation's accuracy the empirical way:
take a problem with a known answer, solve it on a sequence of meshes, and show
the error is already small on a relatively COARSE mesh (and merely polishes as
you refine).  The question this module answers:

> Can that coarse-mesh accuracy be stated directly from the **properties of the
> DtN matrix at various mesh sizes**, instead of from a refinement experiment?

Yes.  The accuracy is a **spectral** fact about the discrete
Dirichlet-to-Neumann operator, visible without ever solving the field problem.

## One operator behind every open-boundary method

Truncate the unbounded exterior at a surface Γ and impose the EXACT transparent
condition there.  That condition is

    ∂u/∂n |_Γ  =  Λ_ext u |_Γ ,

where Λ_ext is the exterior Steklov-Poincaré (Dirichlet-to-Neumann) operator:
it returns, for any boundary trace u|_Γ, the outward normal derivative of the
unique exterior harmonic function that decays at infinity.  Every practical
open-boundary technique realises some operator S_h on Γ targeting this single
Λ_ext -- but in two distinct fidelity classes:

  | Method                | How it realises an S_h on Γ        | Fidelity         |
  |-----------------------|------------------------------------|------------------|
  | Kelvin transformation | FEM on the inverted ball; trace Γ  | →Λ_ext as h→0    |
  | BEM coupling          | S_h = Λ_h = V⁻¹(−½M+K) (bnd only)  | →Λ_ext as h→0    |
  | PML                   | absorbing layer; trace on Γ        | fixed-error model|
  | Infinite elements     | decaying shape functions; trace Γ  | →Λ_ext (basis)   |
  | Asymptotic Robin BC   | S_h = −(1/R) I                     | exact n=0 only   |

Two classes: Kelvin / BEM / (rich-enough) infinite elements are DISCRETIZATIONS
that CONVERGE to the exact Λ_ext for every degree as h→0; PML and the asymptotic
Robin BC are MODEL surrogates with a fixed, mesh-independent modal error (the
plain Robin S_h = −(1/R)I reproduces the n=0 eigenvalue exactly and applies that
same −1/R to ALL higher modes -- exact n=0, wrong n≥1, independent of h).  This
module is about the first class: there the open-BC error is governed by how well
the DISCRETE S_h reproduces Λ_ext **on the modes the solution contains**.

On the truncation surface Γ: the KELVIN transformation IS a spherical inversion
(r ↦ R²/r about a centre), so its Γ is NECESSARILY the inversion sphere of radius
R -- there is no non-spherical Kelvin.  Hence the spherical eigenvalue ladder
−(n+1)/R below is NOT a special case for Kelvin; it is the COMPLETE, general
spectral structure of any Kelvin coupling (the interior physical domain may be any
shape, but it sits inside the Kelvin sphere and the coupling is on that sphere).
BEM / PML / infinite elements, by contrast, admit ANY closed Γ -- on a non-sphere
the relevant spectrum is that surface's own Steklov spectrum, not the −(n+1)/R
ladder.  So for the Kelvin question this module answers, the sphere is the whole
story.

## The CONTINUOUS DtN spectrum is mesh-independent and known in closed form

On a sphere of radius R the spherical harmonics Y_n^m diagonalise Λ_ext:

    Λ_ext Y_n^m = −(n+1)/R · Y_n^m ,     m = −n … n   (multiplicity 2n+1)

    n = 0 monopole   : λ₀ = −1/R
    n = 1 dipole     : λ₁ = −2/R     (validated, sphere_exterior_dtn_eigenvalue)
    n = 2 quadrupole : λ₂ = −3/R
    n = 3 octupole   : λ₃ = −4/R
    …

These eigenvalues are properties of the CONTINUOUS operator -- no mesh, no h.
They are the yardstick: a discrete DtN matrix Λ_h is "good for mode n" exactly
when its corresponding eigenvalue matches −(n+1)/R.

## The matrix property that explains coarse-mesh accuracy

Assemble the dense DtN matrix Λ_h on a sphere at several mesh sizes and read off
its eigenvalues (exterior_dtn_spectrum / dtn_spectrum_vs_mesh).  Three facts
fall out (measured numbers in the NUMERICS topic):

  1. **The low modes are already accurate on the COARSEST mesh.**  On the
     coarsest admissible sphere mesh the n=0,1,2 eigenvalues match −(n+1)/R to
     <~0.25 % (the n=2 quadrupole is the worst at 0.21 %), because the
     eigenfunctions (low spherical harmonics) are SMOOTH: a coarse surface mesh
     interpolates them with tiny error, and the Galerkin eigenvalue (which for an
     analytic eigenfunction superconverges as the SQUARED L2 trace error)
     inherits an even smaller error.  This is the coarse-mesh accuracy itself --
     visible in the matrix, before any solve.

  2. **At EVERY mesh size the error is ordered by degree n** (the spectral
     signature).  rel_err increases SMOOTHLY and monotonically with n -- not a
     cliff.  A surface mesh can represent at most ~√ndof harmonic degrees (the
     angular Nyquist ceiling: (n+1)² harmonics in ndof≈1/h² DOFs; ~18 at
     ndof=336), but the per-degree eigenvalue error degrades gradually as a
     smooth h-power with an n-growing constant FAR below that ceiling -- the
     engineering-accurate band (n≤2 at 0.5 %) sits well under the ~18 it could
     nominally resolve.  So a fixed mesh is "accurate up to some low degree,
     progressively worse above."

  3. **Refinement lowers every mode steeply and WIDENS the accurate band.**
     Each refinement level cuts the per-degree error ~×2.6 here (measured rate
     p≈3.9, an order-1 Galerkin eigenvalue superconvergence ~O(h⁴); slightly
     slower, ~×2.5, by n=4 as the eigenfunction gets less smooth relative to h),
     pushing the resolution limit up so more degrees fall below any fixed
     tolerance.  But the low modes were already below engineering tolerance on
     the coarse mesh -- refining buys you HIGHER modes, not materially better
     low ones.

The honest distinction (not "the low modes are mesh-independent"): the low-mode
error *does* fall under refinement, but it starts from such a small floor on the
coarse mesh that there is nothing of engineering value to gain there.  What
refinement actually adds is bandwidth -- accuracy for higher multipoles.

## Why this matches Kameari's observation

For a COMPACT source inside Γ the exterior field is a multipole series whose
n-th term decays like r^{−(n+1)}.  At a truncation radius a few times the source
size the boundary data is dominated by the lowest harmonics; the high-n content
is small.  Therefore:

  * the OPEN-BC part of the error depends on (S_h − Λ_ext) only through its
    LOW-degree block -- the very block that fact (1) says is already accurate on
    a coarse mesh, and fact (2) says is the most accurate part of the spectrum;
  * Kameari's field-refinement curve is the per-degree convergences of fact (3)
    superposed and weighted by the source's multipole content.  Because that
    content is low-degree dominated and those modes are already accurate coarse,
    the open-BC contribution to the field error is already small on the coarse
    mesh and merely polishes.

This is an EXPLANATION (implication), not a strict equivalence -- it rests on
three premises, all of which the spectrum lets you check:
  (i)   the boundary data is low-degree dominated (truncation a few source-radii
        out, and no sharply-structured source sitting near Γ -- see Applications);
  (ii)  the OPEN-BC error is not the bottleneck: the total field error also
        carries the INTERIOR FEM discretisation error, which is typically the
        larger term on a coarse mesh (e.g. ~5 % L2 from the interior FEM solve
        while the dipole DtN eigenvalue is already 0.07 %).  The spectrum isolates the
        boundary contribution from this interior error -- which a field-only
        refinement study conflates;
  (iii) for the Kelvin transform specifically, that its S_h shares the low-mode
        fidelity -- now MEASURED directly (kelvin_dtn_eigenvalue), see below.

With those premises the refinement experiment and the spectral statement are two
views of the same continuous operator.  The spectral one is stronger: it tells
you the answer is good BEFORE you refine -- and exactly which degrees a given
mesh can be trusted for -- and it isolates the open-BC error from the interior
FEM error, without solving the field problem at all.

## The Kelvin bridge (now MEASURED, not just argued)

The Kelvin inversion r ↦ R²/r' maps the slowly-decaying low-n exterior modes to
BOUNDED solid harmonics on the Kelvin ball.  With the 3D conformal weight
(R/r')^(d−2) = R/r' the material modulation carries (see knowledge:
kelvin_transformation):

    u*_n(r') = (R/r') · u_n(R²/r') = (R/r')·(R²/r')^{−(n+1)} Y_n  ∝  r'^{n} Y_n

so the exterior r^{−(n+1)} Y_n becomes the SOLID HARMONIC r'^{n} Y_n on the ball
-- a DEGREE-n POLYNOMIAL.  Order-p Lagrange FEM represents a degree-n polynomial
EXACTLY iff p ≥ n.  That is the mechanism, and it is measured directly by
`kelvin_dtn_eigenvalue` (a volume-FEM solve on the inverted ball, reading the
effective DtN eigenvalue λ_eff = −1/R − ∫_Ω|∇u*|² / ∮_Γ u*² and comparing to
−(n+1)/R).  Measured on a sphere R=1:

  Kelvin closure effective DtN, rel_err vs −(n+1)/R  (volume FEM):

   mode          inverts to | order 1            order 2     order 3
   --------------------------|---------------------------------------
   dipole    n=1   linear    | 0.5% -> 0.03%      --          --
                             | (maxh 0.6->0.25, converges ~O(h³))
   quadrupole n=2  quadratic | 18%  (maxh 0.4)    0.03%       --
   octupole  n=3   cubic     | 40%  (maxh 0.4)    --          4e-5

Reading it: the DOMINANT dipole inverts to a LINEAR field, so even order-1 FEM
nails its DtN eigenvalue on the coarsest mesh -- the residual is only GEOMETRY
(curved-sphere) error, converging fast; nothing about the dipole needs refining.
The quadrupole inverts to a quadratic and needs order ≥ 2; the octupole a cubic,
order ≥ 3.  So the Kelvin closure is coarse-mesh accurate for exactly the low
modes a compact source radiates -- now a MEASURED fact, matching the BEM Λ_h
table, via the shared continuous Λ_ext.

Two complementary mechanisms for the SAME conclusion:
  * BEM Λ_h (boundary): all low SURFACE harmonics are smooth, so a coarse surface
    mesh resolves them; error grows SMOOTHLY with degree n.
  * Kelvin (volume): mode n inverts to a degree-n polynomial; the closure is
    exact up to FEM order, a SHARP threshold at n = p.
Both put the dominant low modes inside the accurate set on a coarse mesh.

This also explains the failure mode at the BOTTOM of the table: a fixed-potential
or single-eigenvalue truncation (asymptotic Robin BC, S_h = −(1/R) I) gets only
the n=0 mode right and mis-handles the dipole and up -- which is why a plain
truncation needs a far larger air box than a Kelvin / DtN closure for the same
accuracy.
"""

DTN_COARSE_MESH_NUMERICS = r"""
# Theory FIRST: why low modes are accurate, and the sufficient-mesh criterion

(The measured table further down VERIFIES this; the theory comes first.)

## Control parameter (heuristic): facets per wavelength,  q = n h / R

A degree-n spherical harmonic Y_n oscillates with angular wavenumber n (surface
Laplace-Beltrami eigenvalue n(n+1)/R^2), wavelength ~ 2*pi*R/n.  A surface mesh of
element size h resolves it with ~ R/(n h) elements per wavelength, so

    q = n h / R   ~   1 / (facets per wavelength of Y_n)

is a useful RESOLUTION HEURISTIC.  CAVEAT: it is ONLY a heuristic -- the
experiment below shows n and h enter the defect with DIFFERENT powers, so the
defect is NOT a pure function of q.  The angular-Nyquist WALL is q ~ 1
(n_Nyq ~ R/h ~ sqrt(N_surf)); beyond it the mode is unresolved.

## Superconvergence law (theory) + what the experiment CORRECTED

The DtN/Steklov eigenvalue is variational, so for a SMOOTH eigenfunction the
Galerkin eigenvalue error is the SQUARED best finite-element approximation error
of Y_n.  Order-p_s surface elements approximate the smooth Y_n to O(h^{p_s+1}),
so the h-law is

    defect_n  ~  h^{2(p_s+1)}  x  C(n),     C(n) growing with the degree n.

MEASURED (dtn_spectrum_vs_mesh, SurfaceL2 p_s=1, R=1; ndof 336/564/882/1080):

  * h-power CONFIRMED:  defect ~ h^{2(p_s+1)} = O(h^4), rate p ~ 3.9
    (n=2: 2.1e-3 -> 8.0e-4 -> 3.2e-4 -> 2.2e-4 as ndof rises).
  * n-power MEASURED ~ n^2 at fixed mesh (per-degree rel_err grows ~x2.5/degree,
    log-log slope ~ 2.1) -- this CORRECTS a naive "pure q^4" reading: the clean
    part is the h^4 superconvergence, the degree growth is a SEPARATE, milder
    ~n^2 (its exact power is a subtle Steklov question; the datasheet carries the
    constants).

So for p_s=1 the empirical law is  defect_n ~ n^2 (h/R)^4  -- NOT n and h together
as q^4.  Theory fixes the h-power 2(p_s+1) and the sub-Nyquist structure; the
spectrum fixes the n-growth and the constants.

## Sufficient-mesh criterion (from the measured law defect ~ n^2 (h/R)^4)

A surface mesh is sufficient for a source whose highest SIGNIFICANT multipole is
n_src, at tolerance eps, when defect_{n_src} <~ eps:

    (h/R)  <~  eps^{1/4} / sqrt(n_src)
    <=>  N_surf ~ (R/h)^2   >~   n_src / sqrt(eps)     (LINEAR in n_src!).

So the surface DOF grows only ~LINEARLY with the top multipole (x eps^{-1/2}).
You resolve the top significant harmonic to a few facets per wavelength -- you do
NOT "refine to convergence".

MEASURED band vs Nyquist (rel_err < 0.5%): the accurate band sits at a small,
roughly FIXED fraction of the angular-Nyquist ceiling --

    ndof  336 :  n_acc <= 2,  n_Nyq ~ 17   (0.12)
    ndof  564 :  n_acc <= 4,  n_Nyq ~ 23   (0.18)
    ndof  882 :  n_acc <= 5,  n_Nyq ~ 29   (0.17)
    ndof 1080 :  n_acc <= 5,  n_Nyq ~ 32   (0.16)

-- ~0.12-0.18 of Nyquist: you never approach the resolution wall.  A
dipole-dominated source (n_src = 1) is already sufficient on the COARSEST
admissible sphere mesh (low-mode error 2.1e-3 there).

# Measured DtN Spectrum vs Mesh Size  (sphere R=1, BEM SurfaceL2 order=1)

`exterior_dtn_spectrum` assembles Λ_h = V⁻¹(−½M+K) and matches its eigenvalues,
degree by degree, to the analytic ladder λ_n = −(n+1)/R.  `dtn_spectrum_vs_mesh`
runs it across mesh sizes.  Measured on a sphere R=1, BEM SurfaceL2 order=1,
intorder=10 (NGSolve 6.2.26xx).  NB the sphere mesh FLOORS at ndof=336 -- every
maxh ≥ 0.5 gives that same coarsest mesh -- so the three columns below are the
three genuinely distinct refinement levels:

  per-degree relative eigenvalue error  |λ_h,n − λ_n| / |λ_n|

   degree n      | ndof=336    ndof=564    ndof=924      refines →
                 | (maxh 0.5)  (maxh 0.4)  (maxh 0.3)
   --------------|--------------------------------------------------
    0  monopole  | 2.6e-04     9.9e-05     3.6e-05
    1  dipole    | 7.4e-04     2.9e-04     1.0e-04
    2  quadrupole| 2.1e-03     8.0e-04     3.0e-04
    3  octupole  | 5.2e-03     2.1e-03     7.8e-04
    4           | 1.1e-02     4.6e-03     1.8e-03
   --------------|--------------------------------------------------
   accurate band | n ≤ 2       n ≤ 4       n ≤ 4
   (rel_err<0.5%)|

Reading the table:

  * **The dipole/quadrupole are already correct to ~0.07 % / ~0.2 % on the
    COARSEST mesh** (left column, ndof=336).  Nothing was refined to get there.
    This single fact IS the coarse-mesh accuracy -- read straight off Λ_h.
  * **Down each column the error grows monotonically with degree n** (≈ ×2–3 per
    degree) -- the spectral signature: accurate at the bottom, degraded toward
    the resolution limit.  A given mesh is trustworthy only up to some degree.
  * **Across each row the error falls steeply under refinement** (≈ ×2.6 per
    level on average, ~7× over ndof 336→924; measured rate p≈3.9, an order-1
    ~O(h⁴) Galerkin eigenvalue superconvergence -- slightly slower, ~×2.5, by
    n=4 as the eigenfunction gets less smooth relative to h) -- so the low modes
    are NOT mesh-independent; they simply start from a tiny floor.  What
    refinement materially adds is BANDWIDTH: the band widens n≤2 → n≤4.

Interpretation: the physically dominant low multipoles are captured on the
coarsest mesh, so an open-boundary closure built on this DtN is already accurate
there -- which matches Kameari's coarse-mesh result, read off the matrix spectrum
instead of a field-refinement study.  The numbers above are the BEM Λ_h spectrum;
the Kelvin-FEM realisation is measured SEPARATELY (kelvin_dtn_eigenvalue; see the
"Kelvin bridge" section in the overview) and shows the same coarse-mesh accuracy
for the low modes via a complementary, polynomial-order mechanism.  Refining is
for sources that carry high-degree content, which the spectrum flags as the
still-inaccurate high-n band.
"""

DTN_COARSE_MESH_API = r"""
# DtN Spectrum API (bem_integral.py)

## exterior_dtn_spectrum -- one mesh

```python
from radia_mcp.radia_ngsolve.bem_integral import exterior_dtn_spectrum

res = exterior_dtn_spectrum(R=1.0, maxh=0.6, order=1, intorder=10, nmax=4)
res["ndof"]          # boundary unknowns
res["modes"]         # list of per-degree dicts:
#   {"n":0, "multiplicity":1, "lambda_exact":-1.0,
#    "lambda_mean":-0.998, "rel_err":2.1e-3, "spread":0.0}
res["eigenvalues"]   # leading real eigenvalues, descending (closest-to-0 first)
```

The matched `modes` are obtained by taking the leading negative eigenvalues
(sorted toward zero first), bucketing them by the analytic multiplicity 2n+1,
and comparing each bucket mean to −(n+1)/R.

## dtn_spectrum_vs_mesh -- mesh sweep + summary metrics

```python
from radia_mcp.radia_ngsolve.bem_integral import dtn_spectrum_vs_mesh

# NB the sphere mesh FLOORS at maxh~0.5 (R=1): use maxh <= 0.5 for distinct meshes
study = dtn_spectrum_vs_mesh(R=1.0, maxh_list=(0.5, 0.4, 0.3),
                             order=1, intorder=10, nmax=4, band_tol=0.005)
study["per_mesh"]                 # list of exterior_dtn_spectrum() results
study["ndof_list"]                # [336, 564, 924]
study["error_table"]              # error_table[n] = [rel_err at each mesh]
study["coarse_low_mode_max_err"]  # max rel_err over n<=2 on the COARSEST mesh
study["low_mode_max_err"]         # ... over n<=2 across ALL meshes (worst low-mode
                                  #     error in the sweep, set by the coarsest mesh)
study["degree_monotonic"]         # per mesh: rel_err strictly increases with n?
study["degree_growth"]            # per mesh: rel_err(top n)/rel_err(n=0)
study["accurate_band"]            # per mesh: highest n with rel_err < band_tol
```

Key summary numbers (and what they assert):

  * `coarse_low_mode_max_err` -- small (~2e-3): the low modes are already
    accurate on the COARSEST mesh.  The headline coarse-mesh-accuracy number.
  * `degree_monotonic` = [True, …] and `degree_growth` ≫ 1 -- at every mesh the
    error is ordered by degree (the spectral signature: accurate low, poor high).
  * `accurate_band`  -- increases as maxh decreases (e.g. 2 → 4 → 4): refinement
    WIDENS the band of well-resolved modes (Kameari's refinement, mode-resolved).
    It does NOT mean the low modes were inaccurate before -- they were already
    inside the band on the coarsest mesh.

## kelvin_dtn_eigenvalue -- the Kelvin (volume-FEM) companion measurement

```python
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue

r = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.4, order=1, dim=3)
r["lam"]       # measured effective DtN eigenvalue of the Kelvin closure
r["lam_exact"] # -(degree + dim - 2)/R   (3D: -(n+1)/R ; 2D: -n/R)
r["rel_err"]   # dipole order=1: ~1e-3 (linear -> geometry-limited, converges)
               # quadrupole order=1: ~0.18 (needs order>=2); order=2: ~3e-4
```

Solves Laplace on the inverted Kelvin ball/disk with the mode's solid-harmonic
datum and reads λ_eff = offset − ∫|∇u*|² / ∮u*², offset = −(dim−2)/R.  Mode n is
exact iff FEM `order >= n` (the dipole inverts to a LINEAR field -> order-1 coarse
accurate).  Fast (no densification) -- the complement of the BEM
`exterior_dtn_spectrum`.

**`dim=2` (circle inversion)** -- the cross-section case for STATIC APPARATUS /
ROTATING MACHINES.  The Kelvin transform is a spherical inversion, so in 2D the
truncation surface is a CIRCLE and there is no conformal prefactor (offset=0):
the eigenvalue is **−n/R** (not −(n+1)/R), and the order-threshold MECHANISM
(order≥n kills the polynomial error; dipole → linear → order-1 coarse accurate) is
the SAME as 3D -- but with NO conformal weight the realized geometry floor is DEEPER
than 3D (measured ~1e-7…1e-9 at order≥n, vs the 3D sphere's ~5-6 digits).  Verified:
2D dipole −1/R to 0.06 % at order 1 / coarse mesh; quadrupole needs order≥2.

## kelvin_twosphere_shell_dipole -- the lab's REAL Kelvin (validation)

```python
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_twosphere_shell_dipole
r = kelvin_twosphere_shell_dipole(maxh=0.25, order=2)
r["rel_err"]   # ~5.2e-2 shell L2 vs exact dipole (mostly interior FEM)
```

The genuine two-offset-sphere periodic-BC Kelvin (Nagamine convention, material
μ'=(R/r')², GND at the Kelvin centre) solved for the shell dipole -- NOT the
single-ball effective-DtN equivalent.  Reproduces the exact u=R_in²z/r³ to ~5.5 %
at maxh=0.3 order=2, matching the analytic-DtN shell solve -- so the real Kelvin
coupling = the analytic exterior = the single-ball effective DtN.  Confirms the
whole picture on the actual implementation.

## kelvin_vs_exact_open_bc_error -- ISOLATE the Kelvin open-BC error from FEM error

```python
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_vs_exact_open_bc_error

ev = kelvin_vs_exact_open_bc_error(R_inner=0.5, R_outer=1.0, maxh=0.4, order=2,
                                   degree=1, kelvin_maxh=0.6, kelvin_order=1,
                                   include_bem=True)
ev["interior_fem_error"]   # ~5.3e-2  shared FEM error (open BC EXACT here)
ev["kelvin_openbc_error"]  # ~1.2e-3  ISOLATED Kelvin open-BC error (FEM error cancelled)
ev["bem_openbc_error"]     # ~1.7e-3  ISOLATED BEM open-BC error (if include_bem)
```

On ONE shared shell mesh it swaps only the Γ operator -- exact-DtN Robin
(λ=−(n+1)/R, zero open-BC error), Kelvin-DtN Robin (λ=`kelvin_dtn_eigenvalue`),
and the BEM-DtN Schur -- so the interior FEM error CANCELS and
``||u_method − u_exactDtN||`` is the method's PURE open-boundary error.  For a
single-mode shell problem the Kelvin closure acts on Γ as multiplication by its
effective DtN eigenvalue, so Robin-with-λ_Kelvin IS the Kelvin field (no separate
Kelvin mesh / periodic BC).  This is the quantitative evaluation of the Kelvin
method's OWN numerical error, separated from the universal FEM error: ~0.1 %,
~45× below the ~5 % interior FEM error even on a coarse Kelvin mesh -- Kameari's
coarse-mesh accuracy as a SEPARATED ERROR BUDGET rather than a refinement curve.

## kelvin_openbc_error_vs_exterior_mesh -- accuracy vs the EXTERIOR mesh size

```python
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_openbc_error_vs_exterior_mesh

sw = kelvin_openbc_error_vs_exterior_mesh(kelvin_maxh_list=(0.7,0.5,0.35,0.25),
                                          maxh=0.4, order=2, degree=1)
sw["interior_fem_error"]  # FIXED ~5.3e-2 (independent of the exterior mesh)
sw["per_mesh"]            # per exterior mesh: kelvin_ndof, kelvin_openbc_error,
                          #   ratio_fem_over_openbc, kelvin_dtn_rel_err
sw["converges"]           # open-BC error non-increasing as exterior refines
sw["always_below_fem"]    # open-BC error < interior FEM error at EVERY exterior mesh
```

In the Kelvin transform the unbounded exterior IS the Kelvin ball, so the
"exterior-region mesh size" is ``kelvin_maxh``.  This sweeps it with the INTERIOR
mesh held FIXED, so the (fixed) interior FEM error cancels and the reported
open-BC error is purely the EXTERIOR-discretisation contribution.  Measured: the
open-BC error converges ~×4/level (1.2e-3 → 7.5e-5 over kelvin_ndof 58→301) while
the interior FEM error stays at 5.3e-2 and dominates (45×→709×).  So a COARSE
exterior mesh already suffices -- Kameari's exterior-refinement accuracy check,
ISOLATED (refine the exterior, watch the open-BC error converge below the fixed
interior floor).  (The Kelvin ball mesh FLOORS at ndof=58 for maxh ≳ 0.5, the
volume analogue of the sphere-surface floor.)

## Cost / limits

`exterior_dtn_spectrum` DENSIFIES the BEM operators (O(ndof²) memory, O(ndof²)
applies) and runs a dense eigensolve, so keep it to ndof ≲ 1200 (order=1, sphere
maxh ≳ 0.3).  It is a diagnostic / explanatory tool, not a production solver --
the production open BC uses the operator directly (fem_bem_schur) without ever
forming the spectrum.
"""

DTN_COARSE_MESH_APPLICATIONS = r"""
# Using the DtN-spectrum view

## 1. Sizing the air box / truncation radius without a convergence sweep

Decide which multipoles your source carries (a compact magnet/coil is dipole-
dominated; a symmetric one starts at quadrupole).  You only need the DtN to be
accurate up to that degree.  `dtn_spectrum_vs_mesh` tells you the coarsest mesh
whose `accurate_band` already covers it -- no field-level mesh-refinement study
required.  This replaces "refine until the field stops moving" with "resolve the
DtN modes the source excites."

## 2. Choosing an open-BC method

The table at the bottom of the OVERVIEW ranks methods by which DtN modes they
capture.  A plain fixed-potential / asymptotic-Robin truncation keeps only n=0
(or n=0,1) and therefore needs a large air box; Kelvin and BEM DtN capture the
whole low-mode ladder and stay accurate on a compact, coarse mesh.  The spectrum
makes the trade-off quantitative.

## 3. Trusting a coarse Kelvin model

If a Kelvin (or BEM) model looks "too coarse to be right," the spectrum is the
reassurance: as long as the low-degree DtN eigenvalues are on the −(n+1)/R
ladder (and they are, even on coarse meshes), the open boundary is faithful for
the multipoles your source actually radiates.  Do NOT refine reflexively --
refine only if your source carries high-degree content (a sharply structured
source near Γ), which the spectrum will flag as an under-resolved band.

For a hard number, `kelvin_vs_exact_open_bc_error` ISOLATES the Kelvin open-BC
error from the interior FEM error (shared mesh, swap only the Γ operator): the
Kelvin closure's own error is ~0.1 %, ~45× below the ~5 % interior FEM error even
on a coarse Kelvin mesh -- so refining the Kelvin region is wasted effort, the
interior is the bottleneck.  This is the quantitative form of "trust the coarse
Kelvin model": Kameari's accuracy as a separated error budget.

## 4. Diagnosing an open-BC bug

If an open-BC result is wrong, check the DtN spectrum first.  Low-mode
eigenvalues far off −(n+1)/R indicate a real operator/sign/scaling bug (e.g. the
P^TΛP sign, the −½M exterior jump, or a wrong R), NOT a mesh-resolution
problem -- because mesh resolution only POLISHES the low-mode eigenvalues from an
already-tiny coarse-mesh floor (a few ×10⁻³); they sit FAR off the ladder only
when the operator is wrong, not when the mesh is merely coarse.  This separates
"the operator is wrong" from "the mesh is too coarse," which a pure
field-refinement study cannot.

## Relation to the radia-ngsolve DtN hierarchy

  | DtN type      | Approximates           | Coarse-mesh behaviour (measured) |
  |---------------|------------------------|----------------------------------|
  | Asymptotic Robin | exterior, n=0 only  | exact n=0, wrong n≥1             |
  | Kelvin transform | exterior, all n     | order≥n kills poly error; then geom floor |
  | BEM DtN (Λ_h) | exterior, all n        | low-n accurate, error smooth in n|
  | SIBC / GIBC   | conductor interior     | curvature-corrected (separate)   |
  | AGE           | annular air gap        | analytic per harmonic            |

The Kelvin transform and BEM DtN share the SAME continuous operator Λ_ext, and
both are MEASURED here to be coarse-mesh accurate for the low modes -- by
complementary mechanisms: BEM by surface-harmonic smoothness (error grows
smoothly with degree), Kelvin by polynomial representability (mode n exact once
the FEM order reaches n; the dominant dipole inverts to a linear field).
"""


DTN_COARSE_MESH_P_METHOD = r"""
# Kelvin is a p-method, not an h-method (measured) -- and the polyhedron question

The coarse-mesh accuracy is usually told as an h-refinement story (Kameari).
But the polynomial-image mechanism makes the Kelvin closure a **p-method**: mode
n inverts to the degree-n solid harmonic, and order-p Lagrange FEM captures that
polynomial image EXACTLY iff p >= n (in the REFERENCE space; on the curved 3D
sphere the realized accuracy then floors at the geometry error, ~5-6 digits -- see
"The REALIZED floor" below).  So the efficient lever is element ORDER, not mesh
density.

## Theory FIRST: regularity decides h vs p -- so this is really an hp question

Whether p or h wins is set by the SOLUTION REGULARITY (classical approximation
theory), NOT by the method name:
  * ANALYTIC (smooth) solution -> p-refinement converges EXPONENTIALLY (~e^{-b p}),
    h only ALGEBRAICALLY (~h^{p+1} = N^{-(p+1)/d}).  p wins.
  * CORNER-SINGULAR solution (u ~ r^s near a re-entrant corner) -> pure p stalls
    at an algebraic rate (~N^{-2 s}); h-grading toward the corner is needed; hp-FEM
    (geometric mesh grading + linear p increase) converges EXPONENTIALLY.  hp wins.

The KELVIN region is the MOST p-favorable part of any open-boundary problem:
  * its image is a degree-n POLYNOMIAL (solid harmonic) -> p does better than
    exponential, it TERMINATES EXACTLY at p = n;
  * its geometry is a smooth sphere (handled by high-order isoparametric elements);
  * it contains NO source and NO singularity -- infinity becomes the smooth GND.
So the Kelvin transform's gift is that it renders the far-field / infinity part
ANALYTIC, removing exactly the difficulty that would otherwise demand h.  In the
exterior, p is decisively advantageous.

The h-vs-p verdict for the WHOLE machine is therefore decided by the INTERIOR
physical region, not by the open boundary:
  * smooth interior (no re-entrant corners)         -> p throughout;
  * slots / re-entrant corners / PM-iron edges      -> singular -> use hp
    (geometric h-grading toward the singularity, p elsewhere).
"Kelvin is a p-method" is precise for the CLOSURE; the full solver is
hp = (high p in the smooth Kelvin ball + smooth interior) + (h-grading toward
interior singularities).  Pure p is NOT globally optimal when the interior is
singular.

Practical caveats:
  * conditioning grows with p and per-element cost is O(p^{2d}); the useful order
    is a moderate sweet spot (~2-4) -- beyond it conditioning erodes the DOF win
    (the order-4 quadrupole below is already conditioning-limited).
  * if the source sits NEAR Gamma (high-n boundary content) p would have to be
    large -- instead enlarge the air-box R to push significant modes to low n and
    keep p moderate (the (R, p) trade).  h enters the EXTERIOR only through this.
  * "advantageous" depends on the metric: p wins on error-per-DOF; low-order h is
    simpler to implement / adapt.  For the smooth Kelvin sphere both a coarse mesh
    and high-order geometry are trivial, so p's win there is clean.

## MEASURED: the hp claim on the L-shape corner benchmark

2D L-shape, 270-deg re-entrant corner, exact u = r^(2/3) sin(2 theta/3)
(singularity lambda = 2/3); H1-seminorm error vs DOF N:
  * h-version (order 1):  fitted rate alpha_h = 0.357 (-> 1/3 asymptotically;
    finest pairs 0.335).
  * p-version (fixed coarse mesh, order 1..6):  fitted rate alpha_p = 0.661 (~2/3),
    with FLAT pairwise slopes (0.65-0.68, no upward trend) = ALGEBRAIC, not
    exponential.
  * ratio alpha_p / alpha_h = 1.85 (~2).
CONFIRMS: a corner singularity caps BOTH h and p at algebraic rates (p ~ twice h);
neither is exponential -- only hp (geometric grading) is.  So REGULARITY, not the
method name, sets the rate; the Kelvin region (analytic) is where p goes exact.
(hp's own exponential rate is asserted from theory, not run here.)

## p-path vs h-path (measured, kelvin_dtn_eigenvalue, sphere R=1, dim=3)

  quadrupole n=2 (image = quadratic):
    p-path (maxh=0.5 FIXED, raise order):
        order 1  ndof  58  rel_err 9.5e-01
        order 2  ndof 278  rel_err 2.5e-03      <- ONE order step = 374x drop
        order 3  ndof 768  rel_err 1.2e-05
    h-path (order 1 FIXED, refine mesh):
        ndof 58 -> 1075   rel_err 9.5e-01 -> 1.4e-02   (algebraic ~O(h^2), never < 1e-2)

  octupole n=3 (image = cubic):
    p-path: order 1/2/3, ndof 58/278/768   rel_err 2.1 / 0.196 / 8.4e-04
    h-path (order 1):    ndof 58 -> 1075    rel_err 2.1 -> 4.4e-02

Reading: p-refinement is **exact-terminating** (the polynomial image is captured
once p >= n); h-refinement is merely **algebraic** (order 1 can never represent
the degree-n image, only approach it ~O(h^2)).  To hit a target accuracy, p uses
**~20-80x fewer DOF** than h for the higher multipoles (quadrupole ~22x, octupole
~77x).  The dominant dipole (linear image) is already exact at order 1.

## The DOUBLE benefit of raising p (subtle, from the implementation)

`kelvin_dtn_eigenvalue` curves the mesh isoparametrically (Curve(min(order+1,3))).
So raising the order improves BOTH the polynomial field representation AND the
sphere geometry -- which is why the p-path is so steep.  Note the tension:

  * AFFINE (flat) elements represent a physical-space polynomial EXACTLY (order>=n
    is then strictly exact) but mesh the sphere as a FACETED polyhedron (geometry
    error);
  * CURVED (isoparametric) elements fix the geometry but BREAK exact polynomial
    representation (a physical-space polynomial is not a reference-space
    polynomial under a nonlinear map), so "order>=n exact" is only APPROACHED.

p-refinement relaxes both at once; h-refinement only chips at geometry, slowly.

## The REALIZED floor: order>=n is exact in REFERENCE space, geometry-limited in 3D

Be precise about what "exact iff p>=n" means.  It is a statement about the
REFERENCE-space polynomial: order p>=n captures the degree-n solid-harmonic image,
killing the POLYNOMIAL-approximation error.  What then REMAINS is the curved-boundary
(isoparametric) geometry error and the conformal-weight quadrature -- and in 3D that
floor is NOT machine zero.  MEASURED at the threshold order p = n (kelvin_dtn_
eigenvalue, coarse mesh, raising order):

   3D sphere :  n=2 at order 2 -> rel 2.5e-3, then ~1e-5 as order rises further
                n=3 at order 3 -> rel 8.4e-4, then ~3e-5
                => realized accuracy FLOORS at ~5-6 significant digits.
   2D circle :  n=2 at order 2 -> rel 5.7e-5, deepening to ~1e-7
                n=3 at order 3 -> rel 1.2e-3, deepening to ~3.6e-9
                => no conformal weight (offset = -(d-2)/R = 0 in 2D) + higher curve
                   order -> the floor is DEEPER.

This is exactly Kameari's empirical "order 3 on the sphere gives 5-6 digits": that
5-6 digits is the curved-GEOMETRY floor, REACHED once p >= n_src (octupole n=3 ->
order 3) -- it is NOT a limit of the multipole/DtN argument and NOT multipole
"reflection".  So the honest p-method statement is: **order >= n removes the
POLYNOMIAL error and drops you onto the geometry floor; curve the geometry (higher
isoparametric order, or a little h) to go below it.**  For the 2D cross-sections this
committee studies (static apparatus / rotating machines) the floor is deeper still --
the closure is even MORE p-favorable in 2D than the 3D-sphere numbers suggest.
(Scripts: examples/kelvin_transformation/DtN_spectrum/p_vs_h_study.py is the 3D
sweep; the 2D numbers are kelvin_dtn_eigenvalue(dim=2).)

## PROOF the floor IS geometry: hold p>=n + mesh, raise ONLY the Curve order

The decisive evidence that the 3D floor is the curved-sphere GEOMETRY (not the
multipole, not the method, not the polynomial): fix the FE order p>=n AND the
(coarse) mesh, and raise ONLY the isoparametric geometry order k (mesh.Curve(k)).
MEASURED (floor_vs_curve.py, coarse maxh=0.5):

   Curve (geometry) order k :   1 (flat)    2          3
   n=2 (p=3) rel_err        :   1.33e-2     3.8e-4     1.30e-5
   n=3 (p=4) rel_err        :   1.48e-2     4.4e-4     3.42e-5

Raising ONLY the geometry order from k=1 (flat polyhedron, ~1% faceting error) to
k=3 drops the error ~1000x to the 5-6 digit floor -- with the polynomial image and
the FE order p UNTOUCHED.  That isolates the floor as curved-sphere geometry: a flat
truncation is ~1% off, an isoparametric (curved) one reaches 5-6 digits.  (Past k>=3
it plateaus ~1e-5: the residual conformal-weight quadrature / energy-quotient limit.)
Script: examples/kelvin_transformation/DtN_spectrum/floor_vs_curve.py.

## Connection to Cauer Ladder Network (CLN): both are spectral closures

The DtN eigenvalue ladder -(n+1)/R is a SPECTRAL object of the same kind as the
Cauer-Ladder-Network (CLN) eigenmode decomposition: each characterises a closure by
how it transmits the modes a source excites.  CLN folds the INTERIOR response by
circuit order {R_n, L_n}; the Kelvin closure folds the EXTERIOR response by element
order p.  Same idea -- decompose the physics into eigenmodes and resolve only the
modes that matter -- applied to the interior network vs the open boundary.

## Corollary: the exterior VOLUME mesh is irrelevant -- only Gamma matters (p>=n)

A sharp, practically important consequence: with the truncation-surface Gamma fixed,
refining the EXTERIOR (Kelvin-ball) VOLUME mesh does NOTHING for the open-BC accuracy.

  * The open-BC contribution enters the field ONLY through the discrete DtN operator
    Lambda_h on Gamma -- a property of the Gamma TRACE (surface mesh + order) and the
    exterior solution, NOT of the volume tessellation.
  * For mode n the exterior solution is the polynomial P_n = rho^n Y_n.  At order
    p>=n it lies in V_h, so by Galerkin orthogonality the FE solution EQUALS P_n on
    ANY volume mesh.  MEASURED: ||u_h - P_n||_L2 ~ 1.5e-15 (machine zero) at every
    density from 107 to 54,784 tets -- the volume solve is already exact; there is
    nothing for refinement to improve.
  * Then the effective DtN
        lambda = -1/R - integral_Omega|grad P_n|^2 / oint_Gamma P_n^2
    is an integral of a FIXED polynomial over a FIXED domain: it depends on Gamma's
    GEOMETRY alone, the volume mesh enters nowhere.  Every lambda change under
    "exterior refinement" is the boundary mesh re-projecting toward the true sphere
    (a SURFACE effect) -- not the interior density.
  * CONTRAST p<n: the polynomial is not captured (||u_h-P_n|| = 1.36 -> 0.028 under
    refinement), so refining DOES help -- but that signals ORDER deficiency, fixed by
    raising p, not by refining the exterior volume.

PRACTICE: do not waste DOF refining the air-box interior.  Spend DOF on (i) order
p>=n_src and (ii) the Gamma surface mesh/geometry; keep the exterior volume coarse.
(Scripts: examples/kelvin_transformation/DtN_spectrum/kelvin_exterior_mesh.py and
kelvin_exterior_mesh3.py.)

## Polyhedral (flat) truncation: the faceting error SCALES WITH MULTIPOLE DEGREE

Does a FLAT polyhedron truncation (Curve 1) work, or must the sphere be exact?
Measured (effective DtN eigenvalue rel_err, flat Curve-1 vs curved):

  dipole n=1 (order 1, LINEAR image):
    flat   (C1)  maxh 0.6->0.14   rel_err 1.14e-2 -> 8.7e-4   (~O(h^2), small)
    curved (C2)                    5.1e-3  -> 6.6e-5
    flat/curved ratio  2.2 -> 13      => dipole is geometry-ROBUST (polyhedron fine)

  quadrupole n=2 (order 2, QUADRATIC image):
    flat   (C1)  maxh 0.6->0.20   rel_err 1.33e-2 -> 2.86e-3  (~O(h^2), slow)
    curved (C3)                    2.55e-3 -> 7.8e-6
    flat/curved ratio  5 -> 369       => quadrupole is geometry-SENSITIVE

So it is NOT that faceting is globally irrelevant -- its IMPORTANCE GROWS STEEPLY
WITH MULTIPOLE DEGREE.  Why (spectral): the polyhedral faceting is a boundary
perturbation that is HIGH spatial frequency (angular degree ~ R/h, the facet
rate) and SMALL amplitude (chord sagitta ~ h^2/R).  By spherical-harmonic
orthogonality such a perturbation couples a mode n mainly to FAR high modes; the
diagonal (eigenvalue) defect it induces is suppressed for low n and grows with n
(higher Y_n carry more high-frequency content that overlaps the faceting).  PLUS
a benign degree-0 part: an inscribed polyhedron sits uniformly inside the sphere
-> an O(h^2) "effective radius" shift that rescales the WHOLE ladder -(n+1)/R_eff
equally (not mode-mixing; removable by area/radius matching).

CONSEQUENCE for p=1: the only mode captured exactly at order 1 is the dipole
(linear image), which is precisely the geometry-ROBUST mode -- so a flat
polyhedron truncation is fine at p=1.  Raise the order to capture higher
multipoles and you must ALSO curve the surface, because those modes are
geometry-sensitive: **p and geometry-order should rise together.**  (In the
two-sphere periodic Kelvin the interior and Kelvin-image share the SAME discrete
surface -- facet-to-facet identification -- so faceting is not an interface
MISMATCH, only "discrete surface vs ideal sphere", i.e. the benign high-n +
uniform-shift error above.)

## Design rule (the inversion of "refine until converged")

Choose element ORDER p >= (source's highest significant multipole); keep the
exterior (Kelvin) MESH coarse; use the air-box RADIUS R only to LOWER the
required p (boundary multipole content ~ M_n / R^{n+1}, so a larger box pushes
significant modes to lower n).  Optimise the PAIR (R, p) for minimum total DOF --
a third axis beyond "enlarge the box" and "refine the mesh".

## Worked (R, p) optimum: sphere vs faceted source (MEASURED, demo_e_optimal_R.py)

Where the (R, p) optimum SITS depends on the source's own multipole content:

  * a DISK / SPHERE source is a PURE dipole on its truncation circle -> the required
    order p never grows with R, so total DOF ~ (R/a)^2 is MONOTONE: the optimum is
    the SMALLEST admissible R (hug the source; the sphere is the minimal-DOF
    truncation -- exactly Kameari's "order 3 on the sphere is enough" regime).
  * a SQUARE / cuboid source carries higher harmonics, so shrinking R raises the
    required p faster than it saves area.  With the DOF proxy (R/a)^2 * p(R)^2 and
    p(R) = ceil(log eps / log(a/R)), total DOF has an INTERIOR minimum at
    R/a ~ 3 (eps = 1e-4..1e-6) -- NOT the smallest box.

So "sphere => minimal R" and "non-sphere => R/a ~ 3" are the two ends of the same
(R, p) optimisation, with the source shape (its multipole content) deciding which.
(DOF proxy = area x order^2; the (a/R)^4 source decay used here is exact -- see the
square-magnet c_n worked example in the DATASHEET topic.)
"""

DTN_COARSE_MESH_FORMULATION = r"""
# Formulation (Omega vs A) and infinity -- a differential-geometry view

Does the open-boundary / DtN treatment depend on the potential variable (scalar
Omega vs vector A)?  And how is infinity best treated?  The de Rham complex
answers both.

## The exterior field is one object; Omega and A are dual potentials

  Lambda^0  --d(grad)-->  Lambda^1  --d(curl)-->  Lambda^2

  * scalar potential Omega is a 0-form: H = -grad Omega (H a 1-form); enforces
    curl H = 0 STRONGLY, div B = 0 weakly.
  * vector potential A is a 1-form: B = curl A (B a 2-form); enforces div B = 0
    STRONGLY, curl H = 0 weakly.

They are the HODGE-DUAL pair; the constitutive law (Hodge star *, the material)
links H (1-form) and B (2-form).

## "The exterior is between Omega and A" = complementary (dual) bracketing

The two dual formulations BRACKET the true energy (classical complementary
variational principles):

    W_Omega  >=  W_true  >=  W_A

so the exact open-boundary solution lies BETWEEN the scalar and vector results,
and the gap  W_Omega - W_A  is the constitutive (Hodge-star) discretisation error
= a GUARANTEED, computable a-posteriori bound.  Running BOTH dual Kelvin solves on
a coarse mesh therefore CERTIFIES the open-boundary accuracy, not just observes
it.  MEASURED (the bracketing PRINCIPLE, Prager-Synge on Poisson -div grad u = f,
u_ex = sin(pi x) sin(pi y), E_ex = pi^2/4): the primal H1 energy and an
equilibrated RT-flux complementary energy bracket the exact energy
E_primal <= E(f_h) <= E_comp on EVERY mesh, gap ~ O(h^2), equilibration residual
~1e-15, and gap >= the true energy error (GUARANTEED bound).  CAVEAT: the strict
bound needs EXACT equilibration div sigma + f = 0; with mixed RT/L2 this holds for
the projected data f_h = Pi_h f, and a RAW flux for transcendental f can dip below
E_ex by the data-oscillation amount (so bracket the consistent-data energy E(f_h),
which -> E_ex).  A direct two-sphere Omega-vs-A Kelvin bracket is the remaining
TODO.

## The DtN gradient block is formulation-INDEPENDENT

On the truncation sphere Gamma, Hodge-decompose the boundary 1-form trace:

    Lambda^1(S^2) = d Lambda^0  (+)  *d Lambda^0  (+)  H^1,     H^1(S^2) = 0.

The EXACT (gradient) block is exactly the scalar Omega DtN -- the SAME ladder
-(n+1)/R -- in EVERY formulation.  A merely adds the coexact (solenoidal) block,
which in a current-free, simply-connected exterior is non-physical gauge and
removable (A reduces to Omega there).  Hence the coarse-mesh / p-threshold result
is formulation-INDEPENDENT on the physical (gradient) block.

MEASURED (partial): the scalar exterior DtN ladder -(n+1)/R is reproduced
(exterior_dtn_spectrum: n=0..3, dipole -2/R to rel_err 7.4e-4).  A-side check on
the dipole A = m x r/(4 pi r^3) (tangential on Gamma): its normal flux
B.n = curl_Gamma(A_t) matches the n=1 zonal harmonic to rel_L2 3.2e-5, and feeding
that harmonic through the scalar exterior DtN returns eigenvalue -2/R (rel_err
2.5e-4).  So the A-side normal flux IS an eigenvector of the SAME scalar DtN with
the SAME -(n+1)/R eigenvalue -- formulation-independence verified for the dipole.
(Assembling a full independent vector H(curl) exterior DtN operator is a TODO.)

SCOPE of "formulation-independent": it is the continuum OPERATOR / eigenvalue
-(n+1)/R that is shared by Omega and A -- they see the one exterior problem.  The
DISCRETE defects are NOT identical: the A-side 3.2e-5 above is the A-discretisation's
OWN approximation error, generally different from the Omega defect at the same mesh.
So "formulation-independent" qualifies the physics (operator, gradient block), not
the numerics (per-formulation discrete defect).  The dual bracket exploits exactly
this: the two formulations' DIFFERING discrete energies bound the one true energy.

## Material modulation = conformal pullback of the Hodge star (depends on degree)

The Kelvin inversion is the SAME conformal map regardless of variable, but the
material it induces is the pullback of *, which depends on the FORM DEGREE:

  * 0-form Omega -> isotropic SCALAR weight   mu' = (R/rho')^{d-2}
        (3D: R/rho';  2D: 1, i.e. ABSENT -> the -n/R vs -(n+1)/R difference)
  * 1-form A     -> anisotropic TENSOR        mu' = (rho'/R)^2 (I - 2 n n^T)
        (Householder) = the Kelvin pullback Jacobian J = (rho'/R)^2 H
        (see knowledge: kelvin_transformation)

The 2D-vs-3D weight difference is the conformal covariance of the 0-form
Laplacian (conformally INVARIANT in 2D, weight d-2=1 in 3D).

## Infinity = one-point conformal compactification

The "right" treatment of infinity is to compactify (adjoin the point at
infinity); the Kelvin inversion realises this, sending infinity to a regular
INTERIOR point (the GND, rho'=0).  The OPTIMAL variable is the one whose
Hodge-star weight is simplest: in a current-free, simply-connected exterior the
SCALAR Omega is canonical and ISOTROPIC -- the differential-geometry reason the
classic "scalar potential in the exterior air, vector potential inside
iron/conductor" hybrid is optimal.  For certified bounds, use BOTH (the dual
bracket above).

## Edge-element A and the point at infinity: let the DtN ABSORB it

A real practical pain: in the A-formulation with EDGE (Nedelec / Whitney 1-form)
elements, the Kelvin map sends infinity to the BALL CENTRE (rho'=0), which is
awkward -- there is no nodal DOF to pin there (edge DOFs live on edges, not
points), and the pulled-back 1-form material is SINGULAR at the centre
(nu' ~ (R/rho')^2 -> infinity).  For nodal Omega this is trivial (Dirichlet GND
at the centre node); for edge A it is not.

DtN dissolves the problem by REMOVING the point at infinity from the
discretisation:

  1. Do NOT mesh the exterior as a Kelvin ball.  Truncate at Gamma and impose the
     exterior operator Lambda_ext as a (nonlocal) boundary condition on Gamma
     (FEM-BEM coupling).  Then there is NO exterior mesh, NO centre point, NO
     singular material -- the decay-at-infinity condition lives INSIDE Lambda_ext.

  2. The exterior is current-free and simply-connected, so its DtN reduces to the
     SCALAR ladder -(n+1)/R (the gradient block; the A-specific coexact block is
     removable gauge).  Couple the interior edge-A to this scalar exterior DtN
     through the surface relation  B.n = curl_Gamma A_t  on Gamma (the normal flux
     is the surface-curl of the tangential A trace -- a natural edge-DOF
     quantity).  This is an A(interior, edge)–Omega(exterior, DtN) hybrid on Gamma.

  3. What "the value at infinity" BECOMES for A: it is the n=0 (monopole) mode.
     For magnetostatics there is no magnetic monopole, so  oint_Gamma B.n = 0
     identically -- the n=0 mode carries ZERO net flux.  So the infinity condition
     for edge-A is exactly the GLOBAL zero-net-flux constraint on Gamma, which the
     DtN carries in its monopole block and which is a clean CONSTRAINT on edge DOFs
     (a surface-flux condition) -- NOT a point condition at a singular centre.

So for edge-element A, prefer the DtN / BEM boundary coupling on Gamma over the
Kelvin-ball mesh: infinity is integrated out into Lambda_ext, and its only trace
is the (automatically satisfied) zero-net-flux monopole constraint.  The
edge-element infinity difficulty simply disappears.  [Scalar FEM-BEM is
implemented -- see knowledge: fem_bem_schur; the edge-A curl-curl coupling, or the
flux-coupled scalar-exterior variant, is the natural extension and a clean TODO.]

## Unifying statement (FEEC)

The Kelvin pullback is a CHAIN MAP of the de Rham complex (commutes with d up to
the conformal weight), carrying exterior harmonic forms (decaying at infinity) to
POLYNOMIAL (solid-harmonic) forms on the ball.  The order-threshold p >= n is the
statement that the order-p polynomial FEEC de Rham complex reproduces solid
harmonics of degree <= p exactly.  So coarse-mesh / p-method accuracy is a
COHOMOLOGICAL property of the complex -- formulation-independent on the gradient
block, and dual-bracketable for certified error.

## The PERFORMANCE theory is differential geometry too: DtN = sqrt(Delta_Gamma)

The mesh-resolution / superconvergence / h-vs-p analysis (NUMERICS, P_METHOD
topics) was stated in plain approximation theory.  Its differential-geometry
content unifies it with everything above:

  * The exterior DtN is, to leading order, MINUS THE SQUARE ROOT of the boundary
    Hodge-Laplace-Beltrami operator:  Lambda_ext = -sqrt(Delta_Gamma) + lower
    order  (a first-order pseudodifferential operator; principal symbol |xi'|).
    Its eigenFUNCTIONS are the Laplace eigenFORMS on Gamma -- 0-forms Y_n for the
    scalar (Omega) DtN, vector spherical harmonics (1-forms) for the A DtN.  On
    the sphere Delta_Gamma Y_n = n(n+1)/R^2 Y_n, and the exact eigenvalue -(n+1)/R
    equals -sqrt(n(n+1))/R up to the O(1/R) curvature term.  The "facets per
    wavelength" q = n h/R is literally the mesh resolution of the Laplace eigenform
    of eigenvalue ~ (n/R)^2.

  * PERFORMANCE = FEEC eigenvalue approximation of Delta_Gamma.  How well the
    discrete de Rham complex (Whitney / high-order polynomial differential forms,
    Arnold-Falk-Winther FEEC) approximates the Laplace eigenforms governs the
    defect.  The measured superconvergence "eigenvalue error ~ (best
    form-approximation)^2 ~ O(h^{2(p_s+1)})" is exactly the FEEC eigenvalue
    superconvergence, controlled by the bounded COMMUTING COCHAIN PROJECTION pi_h
    (the discrete complex's approximation property) -- not an accident.

  * h vs p = Whitney vs the HIGH-ORDER POLYNOMIAL de Rham complex.  The Kelvin
    image is a POLYNOMIAL FORM (solid harmonic: polynomial 0-form for Omega,
    polynomial 1-form for A); the order-r polynomial complex (P_r^- Lambda^k /
    P_r Lambda^k) CONTAINS it exactly iff r >= n -- so p-refinement is
    exact-terminating in the de Rham sense, Whitney (h) only approximates
    algebraically.  Interior forms of limited regularity (corner singularities)
    need the Whitney complex with geometric grading (hp).  "Kelvin is a p-method"
    = "the exterior image lies in the polynomial de Rham complex".

  * PROBLEM-INDEPENDENCE (the datasheet) is the geometric intrinsicness of the
    Laplace spectrum: Delta_Gamma depends only on the Riemannian geometry of Gamma
    (the sphere -> n(n+1)/R^2), NOT on the source.  The datasheet is the spectral
    approximation property of the discrete de Rham complex for Delta_Gamma -- a
    cohomological / geometric invariant of (Gamma, mesh, order), not of any
    boundary value problem.

So every strand -- spectrum, superconvergence, sufficient mesh, h-vs-p, dual
bracket -- is ONE differential-geometry statement: the discrete de Rham complex
approximating the boundary Hodge Laplacian whose square root is the DtN, with the
Kelvin conformal CHAIN MAP placing the exterior data in the polynomial complex.
"""


DTN_COARSE_MESH_DATASHEET = r"""
# Problem-INDEPENDENT performance: the Kelvin closure has a DATASHEET

The deepest consequence of the spectral view: you can state the open-boundary
method's PERFORMANCE without reference to any particular problem.  A field
refinement study (Kameari) certifies ONE chosen problem; the DtN spectrum
certifies the METHOD.

## The error factorises: method (operator) x problem (source)

For a compact source the exterior field is a multipole series, and on the sphere
the open-BC error is DIAGONAL in that basis.  So the open-boundary contribution
to the error of ANY problem factorises:

    open-BC error  =  SUM_n  c_n(source)  x  defect_n(method)

  * defect_n(method) = the discrete closure's eigenvalue error at degree n,
    |lambda_h,n - (-(n+1)/R)| / |-(n+1)/R| -- a property of the METHOD ALONE
    (order p, mesh h, geometry/curving).  PROBLEM-INDEPENDENT.
  * c_n(source) = the source's n-th multipole content and its r^{-(n+1)} decay to
    the truncation radius.  The ONLY problem-dependent factor.

Kameari's refinement measures the SUM for one problem; the DtN spectrum measures
the per-mode factor defect_n DIRECTLY -- once, for the method -- never solving a
boundary value problem.

## The source factor c_n, worked: a uniformly magnetised square (MEASURED)

c_n(source) is set by the source geometry AND its symmetry -- and symmetry can make
whole BANDS of c_n vanish, which the eigenvalue defect alone can never tell you.
For a 2D uniformly magnetised square (M || x, half-edge a), the edge-charge field on
a concentric circle of radius R contains ONLY harmonics n == 1 (mod 4):

  * the dipole n=1 dominates; n=3 is FORBIDDEN BY SYMMETRY -- c_3 == 0 EXACTLY, not
    merely small.  A single-resolution decomposition shows a spurious tiny non-zero
    n=3 that is pure quadrature noise and DECAYS under refinement (the false-positive
    trap); the true leading correction is n=5.
  * c_5 / c_1  =  (4/15) (a/R)^4   (exact; measured a5/a1 = 1/60 at R = 2a), i.e. the
    problem-dependent source content falls off as (a/R)^4 to the truncation radius --
    NOT the (a/R)^2 one might naively guess.

Lesson: defect_n is problem-INDEPENDENT, but turning it into an actual open-BC error
needs the source's TRUE c_n -- here a quartic onset gated by a 4-fold symmetry.  This
is also why the (R, p) optimum for a faceted source sits at finite R/a ~ 3 (P_METHOD
topic) while a pure-dipole disk wants the smallest R.  (demo_d_multipole_spectrum.py
reproduces the n == 1 mod 4 selection and the 4/15 constant; demo_e_optimal_R.py the
radius optimum -- both in examples/kelvin_transformation/DtN_spectrum/.)

## Why Kelvin gives a clean, ANALYTIC, universal datasheet

Kelvin IS a spherical inversion, so its truncation surface is ALWAYS a sphere;
the exterior DtN is therefore ALWAYS diagonalised by the same spherical harmonics
with the same closed-form ladder -(n+1)/R (3D) / -n/R (2D).  The eigenBASIS and
the exact eigenVALUES are known analytically and are the SAME for every Kelvin
problem -- so the datasheet (per-degree eigenvalue defect vs order/mesh/geometry)
is an absolute specification with a built-in analytic yardstick.  (A BEM closure
on an arbitrary surface is also source-independent but lacks the closed-form
yardstick; PML / asymptotic Robin have fixed modal defects.  Kelvin's spherical
structure is exactly what makes the datasheet clean and universal.)

TWO ERROR SPECIES ON ONE AXIS (important when comparing methods): defect_n means
different things for Kelvin vs PML/Robin.  For KELVIN (and BEM) the closure is EXACT
in the continuum, so defect_n is a DISCRETISATION error that -> 0 as p / h / curve-
order improve (the 5-6-digit 3D floor in P_METHOD is just where geometry takes over).
For PML / asymptotic-Robin the closure is a MODEL approximation, so defect_n has a
non-zero FLOOR that survives h -> 0 (only enlarging the layer / box lowers it).  The
common per-mode axis is still fair for RANKING fidelity, but Kelvin's defect is
order-controllable while PML/Robin's is a modelling floor -- do not read them as the
same kind of number.

## The Kelvin datasheet (measured in this module)

  * eigenvalue-defect spectrum  defect_n vs degree n        -> NUMERICS topic
  * order-threshold             order>=n kills poly error,   -> P_METHOD topic
                                then a geometry floor (~5-6 digits 3D, deeper 2D);
                                (p is the efficient axis, ~20-80x fewer DOF than h)
  * geometry-sensitivity        flat-polyhedron error grows  -> P_METHOD topic
                                steeply with n (dipole robust, quadrupole 369x)
  * formulation-independence    gradient block -(n+1)/R same -> FORMULATION topic
                                for Omega and A
  * scalar readouts             C <- defect_0 (exact), L_ext <- this section below
                                defect_1 (dipole) -- the C/L dual pair

## When the diagonal factorisation BREAKS (the fine print)

The clean diagonal "SUM_n c_n x defect_n" needs four assumptions; state them so the
datasheet is not oversold:
  1. SPHERICAL (circular) truncation -- so harmonics diagonalise Lambda_ext.  On a
     non-spherical air-box the operator is NOT diagonal in any fixed basis; the
     factorisation becomes a (still source-independent) MATRIX, not a per-mode scalar.
  2. LINEAR materials -- superposition of modes.  Saturating iron COUPLES modes; the
     datasheet then applies to the linearised/operating-point problem, not globally.
  3. The source is INTERIOR to the truncation sphere with no near-field overlap
     (a pure multipole region at R); a source hugging Gamma injects high-n content
     that the (a/R)^{n+1} decay no longer suppresses.
  4. The closure error is SEPARABLE from the interior FE discretisation error (the
     error-isolation result).  They are orthogonal axes only to leading order; they
     can couple when the interior mesh is too coarse to feed the closure a clean trace.
Within these (the usual open-boundary magnetostatics setting) the factorisation and
the "certify the method once" claim hold; outside them, "problem-independent" refers
to the CLOSURE operator's modal fidelity, not the realised total error.

## What problem-independence buys you

  1. A-PRIORI certification: certify a (mesh, order, R) BEFORE solving anything --
     no per-problem convergence study.
  2. TRANSFERABILITY: one datasheet covers every problem on that truncation
     (motor, transformer, accelerator magnet); just check the source's multipole
     content against the accurate band -- no re-verification per problem.
  3. METHOD COMPARISON on a common axis: Kelvin / BEM / PML / Robin ranked by the
     SAME per-mode fidelity (the OVERVIEW table), independent of any problem.
  4. SEPARATION OF CONCERNS: method (operator defect)  _|_  problem (source
     multipoles)  _|_  interior FEM error (the error-isolation result).  Three
     orthogonal axes, each measurable alone.

## Computational cost: measured by the DoF INCREMENT (not solve time)

The right cost metric for an open-boundary closure is the DoF INCREMENT it adds to
the FE system -- machine-independent and reproducible -- NOT wall-clock solve time
(solver/hardware/implementation dependent).  Measured that way the Kelvin closure is
CHEAP:

  * The closure adds the inverted Kelvin ball = ΔDoF unknowns.  Meshed as fine as the
    interior it is ΔDoF ≈ N_interior (total ~2x; measured 1.72-2.03x).  But the
    exterior is VOLUME-IRRELEVANT (the Γ-only result), so the ball can be a COARSE,
    Γ-scale sphere.
  * MEASURED (kelvin_openbc_error_vs_exterior_mesh, shell dipole): the COARSEST ball
    -- ΔDoF = 58, about the size of the truncation surface Γ -- already leaves the
    closure error at 1.2e-3, ~1/45 of the interior FE error (5.3e-2).  Refining the
    ball (ΔDoF 58 -> 301) lowers the closure error (1.2e-3 -> 7.5e-5) but it is
    ALREADY a non-bottleneck.  So the open boundary costs ~a Γ-scale coarse ball of
    DoF and NEVER limits accuracy -- the precise sense in which Kelvin is "cheap".
  * The DtN spectrum is the MEASURE: the order-threshold (p>=n_src) + the defect law
    (defect_n ~ n^2 (h/R)^4) give the MINIMUM ΔDoF that still reaches a target
    accuracy a priori -- size the cheapest admissible exterior without a sweep.

KEEP KELVIN SPARSE -- do NOT condense it into the explicit DtN operator.  Condensing
the exterior into the Steklov-Poincare DtN removes the ball's DoF but turns Γ into a
DENSE N_Γ^2 clique: MEASURED nnz 10-20x the sparse extension and growing as N^{4/3}
(vs the sparse extension's N).  Condensation pays off only when N_Γ is small, the
exterior is reused across many solves (factor once), or the dense block is
H-matrix/FMM-compressed.  For a one-off sparse FE solve, keep the Kelvin ball sparse.

KELVIN IS A SPARSE STAND-IN FOR THE BEM EXTERIOR DtN (the dual of "keep Kelvin sparse").
BEM and Kelvin DISCRETIZE THE SAME Λ_ext on the SAME Γ; BEM is exact in the exterior but
its DELIVERED accuracy is capped by the Γ trace it shares with Kelvin -- the exterior
treatment buys nothing the surface cannot resolve.  So on ONE coarse sphere (336 surface
DoF, demo_n) Kelvin at order p>=n MATCHES/beats the dense BEM per-degree (n=1: 3.3e-5 vs
7.4e-4; n=3: 8.4e-4 vs 5.2e-3) while staying ~93% sparse (~20 nnz/row, ~ms) vs BEM's 100%
dense (nnz=ndof^2, ~67 s).  The one structural difference: the dense BEM matrix encodes
the WHOLE ladder up to Γ's Nyquist at once, whereas Kelvin reaches multipole n~p (raise p
to extend).  Hence Kelvin = "the BEM exterior DtN truncated at multipole p, realised as a
sparse volume" -- a faithful, cheaper replacement for a compact (low-multipole) source.
The user's framing is exact: since BEM discretizes the surface anyway, Kelvin-grade accuracy
at the same Γ is all the exterior-exactness can deliver, so meshing/integrating the exterior
"better" is wasted -- use the sparse volume instead.

KELVIN vs H-MATRIX/FMM (is the sparse volume BETTER than compressing the dense BEM?).
Both beat the O(N^2) dense exterior coupling, by OPPOSITE mechanisms (demo_o):
  * H-matrix/FMM: the dense BEM matrix has admissible (well-separated) off-diagonal blocks
    that are NUMERICALLY LOW RANK (smooth Green's function across the gap) -> store rank r
    -> O(N log N). MEASURED: a 50x56 block between two polar caps (gap 0.91) has rel
    singular values 1, 9.5e-2, 8.1e-2,... -> rank 5 @1e-2 (90% saved), 14 @1e-4. But this
    must be BUILT: cluster tree + ACA/analytic low-rank + singular near-field quadrature +
    an inner V^-1 solve + an indefinite FEM-BEM coupling.
  * Kelvin: the SAME geometric far-DoF pair is STRUCTURALLY ZERO in the sparse volume
    Laplacian (MEASURED: 0 nonzeros in that block; ~9-21 nnz/row). Sparse SPD, local element
    integrals, ~ms assembly, composes with the interior FEM as ONE sparse solve. Even the
    asymptotics favour it (optimal multigrid O(N) vs H-matrix O(N log N)), with no ACA/
    inner-solve constants.
VERDICT: for a sphere-able truncation, Kelvin is the simpler & cheaper data-sparse route -- it
is sparse BY CONSTRUCTION, not by compressing an already-dense object. H-matrix/FMM's genuine
edge is narrower than it first appears: an ARBITRARY closed Γ (not just a sphere) and mature
oscillatory-kernel FMM for high-frequency radiation. So "better than H-matrix" is true within
Kelvin's applicability domain (spherical truncation, static/low-frequency), not universally.

CORRECTION (do NOT credit H-matrix with "exterior field access"): the Kelvin ball IS the
compactified exterior, so Kelvin ALSO gives the exterior field at ANY point -- one sets up and
solves on the ball (the prepared, inverted exterior region) and recovers the physical field at
any exterior x (|x|>R) by the INVERSE Kelvin transform u(x) = (R/rho) u'(R^2 x/rho^2) (3D weight
(R/rho)^{d-2}): a single FE point-evaluation at the mapped interior point x' plus the conformal
weight -- no surface integral, no Green function. VERIFIED (demo_p) from just outside Gamma
(rho=1.2) to far (rho=20): dipole rel err 9e-7..2.4e-5; quadrupole 6e-4..2.3e-2 (the larger FAR
relative error is the ball's ABSOLUTE geometry floor ~1e-5 amplified where the mapped value is
tiny -- far points map near the centre, value ~rho^-(n+1)). So exterior-field availability is
NOT a BEM-only advantage; both BEM and Kelvin provide it. The surviving BEM/H-matrix edge is
ONLY the arbitrary (non-spherical) Gamma geometry and oscillatory-kernel machinery.

UNIVERSAL DtN GENERATOR -> FEM-Kelvin IS A SPARSE BEM APPROXIMATION FOR ARBITRARY SHAPES.
The assembled FEM-Kelvin matrix gives the DtN/Steklov operator on ANY surface in the mesh, not
just the truncation Gamma. Condense (Schur complement = energy quotient of the harmonic
extension, since the FEM solve minimises u^T A u so min over interior dofs = u_b^T S u_b with
S = A_bb - A_bi A_ii^{-1} A_ib the DtN) onto a MATERIAL surface (a body's boundary r=R_in<R_out)
and you get the EXTERIOR-plus-open DtN there. VERIFIED (demo_q): the ladder is (n+1)/R_in and is
INDEPENDENT of R_out -- p-converged to ~1e-6 across R_in=0.3/0.4/0.5/0.7, and the Gamma value
(n+1)/R_out (R_in-independent 2,3,4) is excluded everywhere. Adversarially confirmed by a
self-consistent two-region weighted-harmonic solve: the natural weak interface coupling
(duA/dr|_R + (R/s)^2 duB/ds|_R = 0) drives the GROWING-mode coefficient to exactly 0, so the
condensed shell field IS the true infinite-domain decaying harmonic (R_in/r)^(n+1) and R_out
cancels completely -- "infinity baked in by Kelvin; Gamma is just an arbitrary glue radius."
  DERIVATION TRAP (do NOT verify the pullback by naive radial bookkeeping): treating
  mu=(R/s)^2 as a pointwise weight on the spherically-split energy and guessing the ball image
  as s^n gives WRONG residuals (e.g. 2/Rin + Rin^2/R^3 = 4.125 instead of 4.0), and a 3D
  Monte-Carlo is contaminated near the inversion core. The regular weighted-harmonic radial
  mode is s^(n+1) (from div((R/s)^2 grad)=0 -> exponents -n, n+1), and the regions couple by the
  NATURAL weak flux-sum, not a hand-guessed strong flux equality. With both correct -> (n+1)/Rin
  exactly, R-independent, matching the FEM to 1e-6.
The condensation surface need NOT be a sphere (only the OUTER Kelvin truncation must be): a CUBE
body (demo_q) gives a well-posed (SPD Cholesky succeeds) dipole DtN ~5.34, in the O(1/size)
window between its face-radius and corner-radius sphere DtN. THEREFORE FEM-Kelvin on an arbitrary
body surface IS a BEM approximation: it discretizes the SAME exterior-plus-open DtN operator that
BEM assembles densely via the Green function -- here from a SPARSE volume matrix, infinity exact,
no Green function / singular near-field quadrature. Gamma need only be a sphere ENCLOSING the
(arbitrary) body, and the exterior field anywhere is recovered by inverse Kelvin (demo_p). The
ONLY residual BEM-specific niche is oscillatory high-frequency kernels (mature FMM, volume
dispersion) and high-aspect-ratio bodies where a bounding Kelvin sphere wastes air-layer volume.

FEM-KELVIN = THE SPARSE FACTORIZATION OF THE DENSE GREEN'S-FUNCTION KERNEL ("trace the DtN
backward"). The dense BEM DtN kernel Lambda IS the Schur complement of the sparse FEM-Kelvin
matrix A (eliminate A's volume DoFs -> dense Lambda); equivalently A is a SPARSE representation
of the same operator obtained by ADDING auxiliary volume DoFs (thin air layer + Kelvin ball).
The trade is "DoF UP, cost DOWN": MEASURED (demo_r) on one Gamma (R=1, maxh=0.5, both reproduce
the ladder) -- dense BEM 336 DoF, 100% fill, 70 s (singular Green-function quadrature) vs sparse
FEM-Kelvin 768 DoF, 4.7% fill (~36 nnz/row), 0.008 s (local grad-grad): DoF x2.3 UP but assembly
~8800x FASTER and fill 100%->4.7%, because each added volume DoF couples to only ~20-40 neighbours
(O(N) local) instead of all (dense O(N_S^2)). This is a DISTINCT sparsification from H-matrix/FMM:
H-matrix keeps the SAME DoF and low-rank-compresses the dense off-diagonal blocks (data-sparse);
FEM-Kelvin ADDS DoF for STRUCTURAL sparsity (genuine zeros), stays SPD, needs no cluster tree /
ACA / singular quadrature / inner V^-1. Cheapest route to the exterior DtN for a sphere-enclosable
body: keep the sparse FEM-Kelvin matrix and never form Lambda.

LIMIT EQUIVALENCE & SURFACE-RATE-LIMITING. In the refined limit the dense BEM DtN and the
FEM-Kelvin DtN are EQUIVALENT (two discretizations of the SAME continuous operator) and BOTH are
rate-limited by the TRUNCATION-SURFACE resolution -- the FEM-Kelvin volume is Galerkin-exact at
p>=n (the interior is irrelevant; kelvin_exterior_mesh.py) and the BEM kernel quadrature is
converged, so neither the volume nor the quadrature is the bottleneck. MEASURED (demo_s, dipole,
matched order 1, refine the surface): maxh 0.50 -> BEM 7.4e-4 / Kelvin 5.1e-3; maxh 0.40 -> BEM
2.9e-4 / Kelvin 1.4e-3 -- both fall at ~the same rate as the surface refines. The order-1 ~7x gap
is a method CONSTANT (order-1 dual trace space vs the curved-geometry floor), not a rate; the
difference that SCALES is COST (sparse ~ms vs dense ~minutes), not accuracy.

HONEST PRIOR ART / POSITIONING (literature scan 2026-06-14; do NOT overclaim novelty). The pieces
of this DtN/spectral/sparse-factorization view are ESTABLISHED -- but for OTHER closures, not the
Kelvin inversion transform:
  * "All open-BC closures discretize one exact exterior DtN/Steklov-Poincare operator" + the
    -(n+1)/R ladder: established for infinite elements / PML / local-NRBC / BEM (Keller-Givoli JCP
    1989; Givoli reviews 1999/2004; Demkowicz-Ihlenburg 2001; Hohage-Lehrenfeld-Preuss "Learned
    infinite elements" SIAM JSC 2021; Gander Acta Numerica 2022). The Kelvin transform is NOT in
    that list -- the magnetics-Kelvin literature treats it purely as a geometric coordinate map
    (Freeman-Lowther 1989; Stochniol 1992; Brunotte-Meunier 1992).
  * "Dense BEM Lambda = Schur complement of the sparse FE matrix (eliminate interior), vs dense
    BIE": near-verbatim in Demarcke-Rogier (IEEE AWPL 2011); the FE-Schur / BEM Steklov-Poincare
    spectral equivalence underlies BETI (Langer-Steinbach 2003) and Hsiao-Khoromskij-Wendland; the
    magnetics sparse-FE-open-boundary-vs-dense trade is Lowther-Freeman-Forghani 1989. The
    H-matrix/FMM "same-DoF low-rank compression" contrast is drawn by Ying-Engquist "Compressed
    ABCs" 2014. So the operator identity and the sparse-vs-dense trade are NOT new.
  * Material-surface DtN on an arbitrary body (demo_q's claim D): the Scaled Boundary FEM is the
    strong prior art -- Wolf-Song CMAME 2001 (fundamental-solution-less BEM-equivalent DtN on an
    arbitrary boundary), magnetostatic SBFEM by Birk-Reichel-Schroeder CMAME 2022. SBFEM yields a
    DENSE boundary operator by radial scaling; our route keeps the exterior SPARSE and condenses on
    demand -- a variant, not a new idea.
DEFENSIBLE NOVELTY (claim only this, with the hedge): instantiating the operator/spectral/sparse-
factorization view FOR THE KELVIN inversion boundary specifically (hitherto only geometric),
explaining Kameari's empirical "coarse exterior mesh suffices" as a DtN-SPECTRAL property (FE order
= multipole reach), with a measured sparse-vs-dense (fill/time) comparison, for the static-apparatus
community. Required hedge in any write-up: "we make this view explicit for / instantiate it for the
Kelvin transform; we do NOT claim the DtN unification, the eigenvalue ladder, or the FE-Schur=BEM
equivalence as new." Best single evidencing figure: eigenvalue-ladder overlay (analytic -(n+1)/R vs
dense BEM vs FEM-Kelvin Schur) on one surface -- coincide to n~p, then peel off (shows same operator
+ why coarse suffices + Schur=dense-DtN in one plot). Closest single prior work: Demarcke-Rogier 2011.

CORRECTION (2026-06-14 deep-read of the flagged PDFs -- the note ABOVE OVERCLAIMED: "Kelvin
hitherto only geometric" is WRONG). Two closest-competitor papers the first scan missed:
  * Nabizadeh, Ramamoorthi & Chern, "Kelvin transformations for simulations on infinite domains,"
    ACM TOG 40(4):97 (2021): the KELVIN transform IS already used for unbounded FEM (Poisson/
    Laplace/Helmholtz incl. EM) -- so "Kelvin for open-boundary FEM" is NOT new (graphics audience;
    no DtN/spectral/sparse-vs-BEM framing).
  * Boulmezaoud, "Inverted finite elements," ESAIM:M2AN 39(1):109-145 (2005); Boulmezaoud & Kaliche,
    "Stray field computation by inverted finite elements," arXiv:2301.10544 (2023): inverted-FE
    (Kelvin-type) truncation-free MAGNETOSTATICS, and the harmonic (r/r0)^{-(l+1)} decay appears as
    an analytic benchmark (their Appendix A) -- but NOT framed as the DtN spectrum, and with no
    Schur=BEM identity and no sparse-vs-dense table. THIS IS THE CLOSEST COMPETITOR (full text read).
What is now KILLED as individually-new (cite + disclaim): the FE-Schur=dense-DtN/BIE identity
(Knockaert 2008; Demarcke-Rogier 2011; Hohage-Lehrenfeld-Preuss 2021 "Schur complement = best DtN
approximation"); the exterior-DtN/static-stiffness reduction for magnetostatics (SBFEM Wolf-Song
2001; Birk-Reichel-Schroeder 2022 -- but DENSE boundary operator via RADIAL SIMILARITY scaling, NOT
Kelvin inversion); Kelvin/inverted-FE for unbounded fields (Nabizadeh 2021; Boulmezaoud 2005/2023);
the harmonic decay ladder (Boulmezaoud; exterior-Steklov literature).
NARROWED DEFENSIBLE NOVELTY = the CONJUNCTION only: (a) explicitly identifying the Kelvin-inverted
exterior FE block AS a discrete exterior DtN with the -(n+1)/R ladder, (b) a MEASURED sparse-FE-
Kelvin vs dense-BEM fill/time benchmark, (c) for the static-apparatus magnetics-audit audience. No
single source does this conjunction (Boulmezaoud has (a)'s harmonic decay but not the DtN/Schur/
sparse-vs-dense packaging). Frame the paper as a UNIFYING REINTERPRETATION + measurement, NOT a new
method. RESIDUAL RISK: Birk-Reichel-Schroeder 2022 (CMAME 399:115362) full text was paywalled
(abstract only) -- get the PDF before any journal submission to confirm it has no Kelvin link / no
comparable fill-time benchmark. Must-cite (honest): Knockaert 2008; Demarcke-Rogier 2011; Hohage-
Lehrenfeld-Preuss 2021; Wolf-Song 2001; Birk-Reichel-Schroeder 2022; Nabizadeh-Ramamoorthi-Chern
2021; Boulmezaoud 2005 + Boulmezaoud-Kaliche 2023.

NARRATIVE SPINE (how to TELL the contribution): this is an approach to ACCELERATE open-boundary
BEM by FEM-IZING it. The dense open-boundary BEM exterior operator (Steklov-Poincare DtN) is
replaced by a SPARSE Kelvin-inverted volume FE block that computes the SAME operator -> sparse
assembly (no Green function, no singular near-field quadrature) + sparse SPD solve = fast (demo_r:
dense BEM 100% fill / ~70 s -> sparse FEM-Kelvin 4.7% fill / 8 ms on one surface; same DtN spectrum,
fig_dtn_overlay). IMPORTANT honest nuance for this spine: "FEM-izing BEM" is itself an ESTABLISHED
DIRECTION -- SBFEM is literally a "fundamental-solution-less boundary-element method" (Wolf-Song
2001; magnetostatic Birk 2022), FE-Schur = the dense BIE operator (Demarcke-Rogier 2011; Knockaert
2008), BETI (Langer-Steinbach). BUT every one of those FEM-izations still PRODUCES A DENSE boundary
operator (SBFEM radial similarity scaling -> dense, fully-populated; FE-Schur condensation -> dense
N_Gamma^2 clique). The KELVIN twist that survives as the contribution: it stays SPARSE -- never form
the dense boundary operator at all; keep the inverted-exterior FE block sparse and condense only on
demand. So state the contribution as: "a SPARSE, spatial-inversion FEM-ization of the open-boundary
BEM operator that avoids the dense boundary matrix entirely, with a measured fill/time speedup and
the same DtN spectrum, for static apparatus" -- NOT "FEM-izing BEM" in general (that is SBFEM et al.).

EXTERIOR-WITH-MATERIAL = the genuine BEM differentiator (and a scope beyond Nabizadeh/Boulmezaoud).
Because the FEM-ized exterior is a MESHED (inverted) FE region, any exterior material mu(x) is just a
transformed coefficient there -- an INHOMOGENEOUS exterior or a far material body is captured
natively, with NO Green function. The free-space-Green BEM cannot do this (a non-constant exterior
breaks the free-space kernel); the vacuum-exterior Kelvin/inverted-FE papers (Nabizadeh 2021;
Boulmezaoud 2005/2023) do not formulate it either. VERIFIED (demo_t): (1) uniform mu_ext -> body
DtN scales as mu_ext*(n+1)/R_in (mechanism; NB a CONSTANT mu keeps the exterior Laplacian so BEM with
a scaled kernel can also do this trivial case -- do NOT cite uniform mu as the BEM differentiator);
(2) INHOMOGENEOUS exterior -- a magnetic shell mu_s in a real sub-shell, vacuum elsewhere -- shifts
the body-surface DtN off the vacuum ladder and FEM-Kelvin MATCHES the layered-sphere analytic
transfer relation (mu_s=5 -> 5.268 vs 5.268; mu_s=20 -> 6.414 vs 6.414; rel ~5e-5..3e-4). THIS is the
real differentiator: surrounding iron / tank / layered media in the exterior are carried by the
sparse inverted FE, which a free-space-Green BEM cannot represent. PRIORITY NOTE: the authors
(Sugahara et al.) state their own Kelvin/inverted-FE-with-exterior-material work PREDATES Nabizadeh
2021 -- if a dated publication exists it is the authors' OWN prior art (so Nabizadeh is parallel/
later, NOT prior art over them) and the exterior-material formulation is their scope advantage; the
citation (date/venue, likely a Japanese IEEJ venue) must be pinned down before final positioning.

SELLING POINT (sharpened): "kernel approximation WITHOUT the Sommerfeld Green's function." For an
INHOMOGENEOUS/LAYERED exterior a boundary-integral (BEM) treatment needs either the Green's function
OF THAT layered medium -- the Sommerfeld (layered-media) Green's function, whose construction via
Sommerfeld integrals is expensive and a subfield in itself -- OR a volume integral equation over the
inhomogeneity (which revives a dense-coupled volume unknown). FEM-Kelvin needs NEITHER: the layered/
inhomogeneous exterior is just a transformed coefficient in the sparse inverted FE, so it approximates
the exterior BEM KERNEL (operator) without ever forming the hard kernel. demo_t is exactly this proof
(layered magnetic shell exterior matched to the layered analytic with NO Green's function). State the
payoff as: FEM-ization sidesteps the difficult BEM kernel (Sommerfeld/layered Green) in the one regime
where BEM is hardest -- inhomogeneous exteriors -- while staying sparse. RELATED-WORK citing: even with
the authors' priority, cite Nabizadeh 2021 / Boulmezaoud 2005-2023 / SBFEM as related/parallel work
(noting priority), NOT as preempting prior art.

CONSTRUCTIVE COROLLARY (demo_u, verified 2026-06-14): FEM-Kelvin not only AVOIDS the layered Green's
function -- it CHEAPLY GENERATES its ACTION. The sparse SPD Kelvin-FEM matrix A IS the discrete
Green's OPERATOR of the (stratified) exterior: A^{-1} applied to a source = the field everywhere,
including the far field read back through the inverse Kelvin map. So the field of a source radiating
THROUGH a magnetic shell -- what a BEM would integrate the (expensive) Sommerfeld/layered-media
Green's function to obtain -- is ONE sparse solve (~43 nnz/row, fill < 0.1%): no Sommerfeld integral,
no special functions, no branch cuts / surface-wave poles, no singular quadrature; the stratification
rides as an FE coefficient and infinity is the Kelvin centre. VERIFIED: a Dirichlet source Y_n on a
body r=a, shell mu_s in [0.7,0.9] vacuum elsewhere, reproduces the analytic layered radial profile
R_n(r)=A r^n + B r^-(n+1) at ALL radii -- through the shell AND out to r=10 via inverse Kelvin -- to
rel ~1e-4..3e-3 (n=1,2; mu_s=1,10,50). HONEST BOUNDS (state these, do NOT overclaim): (1) cheap = the
OPERATOR/ACTION (one solve) and any single value; the FULL dense kernel (all source-observer pairs) is
still O(N^2) entries == the Schur complement of THIS matrix (demo_r) -- Kelvin-FEM makes each entry
cheap, NOT the object smaller. (2) it is a NUMERICAL operator: accuracy = FE order x surface reach (the
-(n+1)/R peel-off datasheet), not a closed-form special function. (3) the genuinely HARD Sommerfeld
case is PLANAR stratification (half-space): under inversion the flat interfaces become spheres through
the inversion centre -- representable but geometrically heavier (Sugahara's extended-Kelvin / radiating
regime, sugahara2025). (4) the inverse map is CONVENTION-dependent: the PERIODIC-GLUE construction
(potential continuous across Gamma) reads u_phys(x)=u'(x') with NO (R/rho) weight, whereas a standalone
compactified ball storing the *transformed* field uses u_phys=(R/rho)u' (demo_p); and a free single
layer with only a one-vertex gauge blows the constant mode up ~1e5x, so inject the source as a
Dirichlet body surface. FRAMING for the audience: "FEM-Kelvin is a sparse numerical GENERATOR of the
layered (Sommerfeld-type) Green's operator" is the positive way to state demo_t's BEM differentiator.

DIRECT ASSEMBLY OF THE MATERIAL-LOADED EXTERIOR DtN MATRIX (demo_v, verified 2026-06-14). "Can the
matrix be made directly?" YES. The dense exterior DtN / Green matrix Lambda on a body surface Gamma is
the SCHUR COMPLEMENT of the sparse SPD Kelvin-FEM operator A: keep the Gamma DoFs, eliminate the rest,
  Lambda = A_GG - A_GR (A_RR)^{-1} A_RG    (ONE sparse factorization of A_RR + a back-sub per Gamma DoF).
Lambda is small, DENSE, symmetric, SPD (on n>=1), size = #Gamma surface DoFs -- a reusable boundary
operator (impedance / macro-element / numerical Green matrix) mapping Dirichlet data on Gamma to the
exterior+open Neumann flux, infinity baked in (Kelvin), NO free-space/Sommerfeld Green's function ever
formed. Because A carries the exterior mu(x) as an FE coefficient, Lambda is the MATERIAL-LOADED
exterior DtN -- the layered/Sommerfeld Green operator AS A MATRIX -- which a free-space-Green BEM cannot
assemble without the (expensive) layered Green's function. VERIFIED: magnetic shell mu_s in [0.7,0.9]
outside a body r=a=0.5 -> Lambda is 452x452 dense (from 12119 sparse DoF, A fill 0.35%, built ~8s,
sym 0.0); its generalized spectrum (Lambda, M_Gamma) reproduces the analytic vacuum ladder (n+1)/a AND
its layered shift, with correct 2n+1 multiplicities:
  mu_s=1 : n=1 4.0005/4.0,   n=2 6.003/6.0,    n=3 8.015/8.0
  mu_s=10: n=1 5.8655/5.8649, n=2 7.608/7.604,  n=3 9.140/9.123
  mu_s=50: n=1 6.937/6.936,  n=2 8.115/8.111,  n=3 9.405/9.387   (rel 1e-4..2e-3 = geometry/order floor).
HONEST SCOPE: n>=1 only (dipole/inductance, quadrupole -- the magnetostatically relevant modes). The
n=0 MONOPOLE is spurious (~0): a single ground POINT has ZERO capacity in 3D H1, so a Gamma-constant
extends as a near-free global constant (energy ~0) instead of the decaying (a/r) monopole; ad-hoc
constant-deflation lands ~3.0, not 2/a. Net-charge/monopole problems need a FINITE ground (small grounded
ball at the Kelvin centre) or the demo_p weighted convention; magnetostatic apparatus has no monopole,
so n>=1 is the whole spectrum. The full dense kernel is still O(N^2) entries (== this Lambda) -- cheap to
GENERATE (sparse factor + back-subs, no Green fn / singular quadrature), NOT made smaller.
PAPER FRAMING (constructive headline): "a directly-assembled, MATERIAL-AWARE exterior DtN (Green) matrix
from sparse Kelvin-FEM, with a closed-form DtN-spectrum datasheet for its rank/accuracy, valid for
inhomogeneous exteriors where the free-space-Green BEM fails." Schur=DtN itself is standard (Demarcke
2011 / Knockaert 2008 / SBFEM); the genuinely-new COMBINATION = exterior-MATERIAL loading (sugahara2022)
+ the -(n+1)/R spectral a-priori + the sparse-Kelvin realization. Cite those as related work, claim the
combination, not the bare Schur=DtN.

NON-SPHERICAL BODY: DIRECT DtN MATRIX VERIFIED BY SYMMETRY (demo_w, 2026-06-14). The BEM-competitive
claim is ARBITRARY geometry. The SAME Schur condensation assembles Lambda for a CUBE body. No analytic
ladder exists for a cube, so verify by the body's OCTAHEDRAL group O_h instead of by a value:
  l=1 -> T1u                  (stays 3-fold)        = the dipole
  l=2 -> E_g(2) + T2g(3)      (sphere 5-fold SPLITS 5->2+3) = the quadrupole
  l=3 -> A2u(1) + T1u(3) + T2u(3).
MEASURED (cube half-side 0.4; Lambda 389x389 dense, sym 0, built ~1s): dipole 3.357 x3; quadrupole
4.531 x2 (E_g) + 5.730 x3 (T2g) -> the 2+3 split is an analytic-value-FREE proof that Lambda is the true
exterior DtN (it carries the body's symmetry). CALIBRATION: same harness on a SPHERE recovers (n+1)/a,
mult 2n+1. p-CONVERGENCE: cube dipole 3.545 (p2) -> 3.365 (p3) -> 3.342 (p4) = a well-defined operator,
not a mesh artefact. MATERIAL: an exterior shell mu_s in [0.75,0.9] shifts the cube dipole 3.357 ->
4.370 (mu_s=10) -> 5.007 (mu_s=50). (Tet-mesh symmetry is imperfect so degeneracies are approximate
~1-2%; a symmetric mesh would sharpen them.) KEY POSITIONING (refines demo_o's hedge): the BODY surface
Gamma may be ARBITRARY while the Kelvin INVERSION surface stays a free spherical truncation around it --
so the one advantage BEM/H-matrix retained over a sphere-able Kelvin truncation (arbitrary Gamma,
demo_o) is MET HERE by condensing onto the arbitrary body surface; the exterior stays a sphere-able
Kelvin ball. Net: a sparse, material-aware, BEM-STYLE boundary operator for an arbitrary body, no
Green's function ever formed -- BEM's arbitrary-Gamma parity PLUS exterior material BEM cannot carry.

SOMMERFELD HALF-SPACE ISOMORPHISM (demo_x, verified 2026-06-14). The GENUINE Sommerfeld case is a
PLANAR half-space (not the spherically-layered demo_t/u). Put the Kelvin inversion CENTRE ON the
interface (origin on z=0): a plane through the centre maps to ITSELF; Kelvin inversion x->R^2 x/|x|^2
PRESERVES sign(z) so the two media z>0(c1)/z<0(c2) are preserved; infinity->the centre. Hence the
UNBOUNDED half-space problem is CONFORMALLY ISOMORPHIC to a BOUNDED two-region ball problem with the
SAME planar interface and media; the Sommerfeld radiation/decay condition becomes a single point
condition at the centre. The Kelvin-FEM IS the Sommerfeld problem in inverted coordinates -- isomorphic,
not merely analogous. VERIFIED vs the closed-form STATIC Sommerfeld (image) Green's function: a
net-neutral vertical dipole (+q@0.40,-q@0.20) in medium 1 (c1=1) over c2=4 -> image charges
beta=(c1-c2)/(c1+c2)=-0.6 in z>0, transmitted 2/(c1+c2) in z<0; Kelvin-FEM matches in BOTH media AND
the far field (inverse Kelvin) to max rel 2.6e-3 (order 3, IfPos interface, maxh 0.13).
CRUCIAL GOTCHAS (each a ~1e7x blow-up if missed): (1) the open Kelvin compactification has a CONSTANT
near-null mode (the non-decaying gauge); a single ground POINT has ZERO capacity in 3D H1 and cannot
hold it -> use a mean-zero NumberSpace constraint (int_Omega u dx=0; gauge offset then ~1e-3, well
conditioned, umfpack) OR a finite ground ball; a bare point ground leaves the gauge ~1e4x the signal
(usable only if the source is exactly neutral, then subtract the offset). (2) the SOURCE must be EXACTLY
discretely neutral -- normalise by the MESHED ball volumes (net charge ~1e-15), else a tiny net-charge
residual excites the n=0 mode ~1e7x. (3) keep the truncation sphere as ONE face (single periodic
Identify) and carry the interface by IfPos(z,c1,c2); a geometric z=0 hemisphere split = two Identify
calls glued imperfectly -> a weakly-coupled Kelvin ball = yet another near-null mode that blows up.
HONEST: this is the STATIC Sommerfeld (image limit). The full WAVE Sommerfeld integral (branch cuts,
surface-wave poles) is the time-harmonic case = the authors' extended-Kelvin radiating regime
(sugahara2025, Maxwellian PML in the Kelvin-mapped exterior); the centre-on-interface isomorphism is
IDENTICAL, only the per-region operator becomes Helmholtz. PAPER VALUE: the strongest BEM-differentiator
statement -- Kelvin-FEM is provably the SAME operator as the Sommerfeld layered-media problem
(isomorphic), realised sparsely and WITHOUT the layered Green's function that is the hard part of BEM
for stratified media.

SOMMERFELD KERNEL = the open-math BEM reference (demo_y, verified 2026-06-14; NOT in NGSolve/ngbem,
which has free-space Laplace/Helmholtz kernels only, nor in the Radia core). A static layered-media
Green's function for a planar stack, source/obs in the top half-space:
  G = 1/(4 pi c0) [ 1/|x-x'| + INT_0^inf R(k) e^{-k(z+z')} J0(k rho) dk ]
-- the reflected part is the STATIC SOMMERFELD INTEGRAL, R(k) the layered reflection coefficient from
the interface recursion R_i=(r_i + R_{i+1}e^{-2k t_{i+1}})/(1 + r_i R_{i+1}e^{-2k t_{i+1}}). For two
half-spaces R=r01 const -> the single image (=demo_x); for a SLAB R(k) is k-dependent -> a true
Sommerfeld integral = an infinite image series. VERIFIED: numerical integral == closed image series to
~1e-16 (slab), and r12=0 / t->0 limits collapse to the single two-media image exactly. WHY for the
paper: it is the independent, open REFERENCE the Kelvin-FEM (demo_x) is isomorphic to -- the one a
single image cannot give for >2 media -- so it certifies the isomorphism for genuine multilayer stacks.
PURPOSE (state clearly, avoid the irony trap): a Sommerfeld kernel is NOT for solving open boundaries
(Kelvin-FEM already does that, AVOIDING the Sommerfeld integral); it is (a) the gold-standard BENCHMARK
for the isomorphism, and (b) a complementary BEM capability for stratified media. GOTCHA: the k->0
spectral boundary layer (width ~1/(2t)) is the classic delicate part of Sommerfeld-integral numerics --
a thick slab's DC limit approaches the 2-media image only as O(1/t) because R(k=0)=r02 for ANY t;
under-resolving k~0 spuriously 'passes' (this is why DCIM / tail-extraction methods exist: Aksun;
Michalski-Mosig). The full WAVE Sommerfeld integral (branch cuts, surface-wave/Zenneck poles) is the
time-harmonic extension. LIT (W:\03..\11_BEM_..\10_sommerfeld_layered): Koh-Yook 2006 (impedance plane),
Chew lectures (Weyl/spectral) [pre-existing]; ADDED 2026-06-14: Sautbekov 2020 (short-wave asymptotic,
arXiv 2005.06857), Lai-Greengard-O'Neil 2015 (hybrid integral representation, arXiv 1507.03491),
Burke-Poggio NEC-2 Part I theory (the classic PUBLIC-DOMAIN Sommerfeld-ground MoM = the user's "NES"/NEC
memory; SOMNEC/SOMNTX). HOME: pure numpy/scipy in radia-ngsolve examples (open math = public-boundary
OK); NOT the C++ Radia magnetostatics core unless a strong reason emerges.

SOMMERFELD LAYERED-MEDIA FORMULATION & NUMERICS -- literature precis (2026-06-14 deep-read of NEC-2
Part I theory [Burke-Poggio], Sautbekov 2020, Koh-Yook 2006, Lai-Greengard-O'Neil 2015, Chew lectures).
Background for deciding whether/how to build a Sommerfeld capability; also sharpens WHY the static
Kelvin-FEM isomorphism (demo_x/y) is the BENIGN limit.
SPECTRAL FORM (Chew): exp(ik0 r)/r = i INT_0^inf (krho/kz) J0(krho rho) exp(i kz|z|) dkrho, kz=
(k0^2-krho^2)^(1/2), physical sheet Im kz>0 (Sommerfeld identity; 2D plane-wave version = Weyl). LAYERED:
split source into TE/TM spectra, insert the generalized reflection coefficient R~^{TE/TM}(krho) (per-
interface recursion) into the integrand. STATIC k0->0: kz->i krho, exp(i kz|z|)->exp(-krho|z|), kernel
-> Lipschitz-Hankel 1/sqrt(rho^2+z^2)=INT e^{-krho|z|}J0(krho rho)dkrho, R~ -> image factor
(eps1-eps2)/(eps1+eps2) = EXACTLY demo_y.
NEC-2 (the user's "NES"/NEC memory; CAVEAT the PDF is a SCANNED image, eqs OCR-reconstructed -- but the
two integrals match the standard half-space forms): fields = free-space + perfect-image CLOSED forms
PLUS two numerical Sommerfeld integrals (J0 kernel)
  U22 = INT_0^inf [e^{-gam2(z+z')}/(gam1+gam2)]            J0(lam rho) lam dlam     (TE-like)
  V22 = INT_0^inf [e^{-gam2(z+z')}/(k1^2 gam2 + k2^2 gam1)] J0(lam rho) lam dlam     (TM-like)
gam_j=(lam^2-k_j^2)^(1/2). Three ground models by cost/fidelity: perfect-image (exact, free); Fresnel
reflection-coefficient image (Rv,Rh; correct only for plane-wave illumination; ~2x); rigorous
Sommerfeld/Norton (U22/V22 for interaction<1 wavelength, Norton asymptotics beyond; wires only; both
points air-side). SOMINT/SOMNEC = PRECOMPUTE the smoothed integrals on a grid in (R1, theta=
atan((z+z')/rho)) after PEELING exp(-jkR1)/R1 and the sin/cos(phi) factors, then bivariate-INTERPOLATE
at fill time (per eval ~0.06s -> one table ~15s). Reusable structure: peel the closed-form free-space/
singular part, tabulate ONLY the smooth ground remainder, interpolate.
THE FOUR NUMERICAL DIFFICULTIES + fixes:
 (1) BRANCH POINTS krho=+-k_j (kz two-sheeted) + log branch at krho=0 (Hankel). Fix: Sommerfeld
     Integration Path (SIP) -- real axis but passing ABOVE -k0, BELOW +k0; vertical cuts (NEC); Im kz>0;
     deform into complex lam to ride exp(-gam2(z+z')) decay; guard the virtual pole crowding +k2 as
     k1->k2.
 (2) SURFACE-WAVE / ZENNECK POLES of R~ near the SIP -> detour (+residue if crossed). SUBTLETY
     (Sautbekov, full-wave): the lossy half-space "surface wave" is a PSEUDO-wave -- the steepest-descent
     contour never captures the R_par pole, so no independent Zenneck residue; it is a near-grazing
     transition (Sommerfeld numerical distance p=k0 rho delta^2<1) handled uniformly by erf / parabolic-
     cylinder (Weber) etalon integrals; a lateral/head wave (total-reflection branch point) decays 1/rho^2.
 (3) SLOW OSCILLATORY TAIL: J0~ (krho rho)^(-1/2)cos -> only conditionally convergent. NEC: Romberg on
     ~0.2pi/max(rho,z+z') sub-intervals + Shanks transform; modern standard = Mosig weighted-averages /
     partition-extrapolation over J0 half-period zeros + analytic tail/quasi-static extraction. Switch
     kernel by geometry: Bessel J0 when rho<=(z+z')/2, Hankel H0^(2) (exp decay) otherwise (incl. z=z'=0).
 (4) NEAR-INTERFACE / small-rho near-singularity + per-pair cost. Modern closed-form route = DCIM
     (discrete complex image method; NOT in the 3 read papers = external knowledge): extract quasi-static
     + surface-wave terms, fit the remaining spectral kernel by a few COMPLEX EXPONENTIALS in kz (GPOF /
     matrix-pencil), apply the Sommerfeld identity term-by-term -> CLOSED-FORM spatial complex images
     (sum exp(ik0 r_i)/r_i, complex r_i). Lai-Greengard-O'Neil alternative: SPLIT G = windowed near-
     singular LOCAL part (physical-space high-order QBX/Nystrom layer potentials) + smooth remainder via
     ONE finite-range Sommerfeld correction whose spectral density decays superalgebraically INDEPENDENT
     of source height (err 1e-10..1e-15). That "physical-space-local + spectral-smooth" split is the SAME
     instinct as Kelvin-FEM (mesh the near/material part, no Green's function) -- a paper framing.
VERIFICATION REFERENCES: Koh-Yook exact closed forms (impedance plane), eta->0 (PEC)/eta->inf limits =
-E1(...) exponential-integral closed forms (clean checks); geometric-optics image = leading asymptotic;
demo_y slab image-series == static Sommerfeld integral to 1e-16.
DECISION CONCLUSION: ALL the hard machinery (SIP, branch cuts, Zenneck/pseudo poles, slow tail, DCIM) is
WAVE-REGIME (k0!=0). The STATIC limit (the SA paper's regime) k0->0 KILLS every one: kz->i krho single-
valued on krho>0 (no branch points), no propagating poles, EXPONENTIAL (not oscillatory) tail -> the
static layered Sommerfeld integral is BENIGN and demo_y already nails it to 1e-16. So a STATIC reference
costs ~nothing (done); a WAVE Sommerfeld capability is a real subproject (SIP + pole tracking + tail
acceleration + DCIM) worth it ONLY for the radiating/extended-Kelvin paper. PAPER POINT: the century of
Sommerfeld-integral difficulty is a WAVE phenomenon; the static stratified problem the Kelvin-FEM is
isomorphic to is its benign limit -- precisely why a sparse real-SPD FEM-Kelvin volume solve can stand
in for it so cleanly. HOME confirmed: Python/numpy/scipy (open math); NOT the C++ Radia core.

FREQUENCY REGIME & THE TWO-PAPER FRAMING (demo_z, verified 2026-06-14). "Is Sommerfeld high-frequency
only? usable at low frequency?" Resolve by separating the KERNEL from the FORMULATION:
 * KERNEL (the layered Green's function) is well-defined at ALL frequencies and NUMERICALLY EASIEST at
   low frequency. VERIFIED: the frequency-dependent half-space Sommerfeld integral vs the PEC exact
   image exp(i k0 R_im)/R_im (R=+1 const -> Sommerfeld identity, exact at any k0): rel err 9.7e-9
   (quasi-static) growing to 3.6e-6 at k0*(z+z')=30 (deep wave) -- the error GROWS with frequency
   because the oscillatory tail IS the difficulty; k0->0 reproduces demo_y (1/R_im, 1.5e-6); a lossy
   half-space (eps_r + i sigma/we) converges at every frequency. So Sommerfeld is NOT inherently high-
   frequency; the 'hard' reputation = the high-frequency oscillatory tail + branch/pole on the SIP.
   Low-frequency layered Sommerfeld is a MATURE field: geophysical EM (CSEM, magnetotellurics), eddy-
   current testing, induction logging -- all quasi-static/diffusive (sugahara2022 is itself eddy-current
   Kelvin).
 * FORMULATION: the famous 'low-frequency breakdown' is an EFIE/MFIE conditioning issue (vector vs
   scalar potential scale as omega vs 1/omega) cured by loop-tree / Calderon -- a property of the
   integral EQUATION, NOT the kernel; and it is ABSENT for the static/quasi-static (Laplace/diffusion)
   formulations the SA paper uses.
THREE REGIMES: (a) DC-static (dielectric/permeable layers): pure-conformal Kelvin isomorphism (demo_x),
benign (demo_y). (b) QUASI-STATIC / eddy-current (low freq, CONDUCTIVE layers, k^2=i w mu sigma): the
SA-relevant regime (transformer tank/core eddy loss; = sugahara2022) -- mesh the conductive layer with
the complex coeff and keep Kelvin for the non-conducting exterior (which stays HARMONIC, so Kelvin still
applies), i.e. our demo_x/v/w machinery extends here directly. (c) FULL-WAVE/radiating (k0 real): hard
Sommerfeld (SIP, Zenneck/pseudo poles, slow tail, DCIM); the Kelvin side becomes extended-Kelvin +
Maxwellian PML (sugahara2025), NOT a pure conformal isomorphism (honest: a truncation-EQUIVALENCE,
weaker than the static operator isomorphism).
TWO-PAPER FRAMING (user, 2026-06-14): the unifying thesis = ACCELERATING / lightening the Sommerfeld-
type LAYERED Green's BEM by replacing the dense special-function integral with a sparse transformed-FE
(Kelvin) solve -- SPLIT BY FREQUENCY because the difficulty (hence the acceleration payoff) differs:
  PAPER L (low-freq / quasi-static = the SA paper, near-complete): Kelvin-FEM as a sparse real-SPD,
    material-aware stand-in for the (benign) quasi-static layered Sommerfeld operator; isomorphism +
    DtN-spectral datasheet; demo_x/v/w/y/z are the evidence; connects sugahara2022.
  PAPER H (high-freq / radiating, a genuine subproject): accelerate the HARD wave Sommerfeld (avoid
    SIP/DCIM entirely) via extended-Kelvin + PML; benchmark against a full wave Sommerfeld reference
    (demo_z is the seed). Caveat: wave = conformal+PML = truncation-equivalence, not a pure isomorphism;
    needs BOTH a wave Kelvin-FEM and a wave Sommerfeld reference built.
NOVELTY DISCIPLINE: 'lighten open-boundary BEM with transformed-FE' is Remacle ~1995 / Lowther 1989 but
for FREE-SPACE exteriors. The defensible NEW angle for both papers = the LAYERED/Sommerfeld-media target
(+ the spectral DtN datasheet, + the wave/PML extension). Re-search the layered-specific literature
(Michalski-Mosig; geophysical-EM transformed-FE) before claiming, but it is a stronger position than the
free-space lightening that Remacle already published.
VERIFIED low-freq (demo_aa): the QUASI-STATIC EDDY-CURRENT (diffusive) conducting-half-space Sommerfeld
integral [R=(krho-u1)/(krho+u1), u1=sqrt(krho^2+i*omega*mu*sigma)] converges and is benign at EVERY
induction number N=sqrt(omega mu sigma)*h (monotone exp tail, no real poles/branch cuts); hits the exact
magnetostatic limits R->0 (N->0, transparent) and R->-1 (N->inf, diamagnetic mirror = -1/R_im, 4.6e-3);
and at moderate-to-high induction (N>=3) reduces to the textbook BANNISTER/WAIT COMPLEX IMAGE (image at
complex depth h+2/gamma, gamma=sqrt(i omega mu sigma); 3e-4 at N=10) -- the low-frequency ancestor of
DCIM, a leading-order approx that degrades at low N (the integral is exact at all N). => the definitive
"Sommerfeld works at low frequency" evidence in the SA-relevant eddy-current regime (sugahara2022).

NON-LAYERED EXTERIOR = the BROADEST differentiator (demo_bb, verified 2026-06-15). The Sommerfeld/
layered Green's function REQUIRES stratified (planar/cylindrical/spherical) symmetry -- that translational
symmetry is what enables the spectral (Hankel) decomposition. For an ARBITRARY inhomogeneous exterior (a
localized inclusion, a curved interface, any mu(x)) there is NO layered Green's function; BEM falls back
to a VOLUME integral equation (the dense volume unknown revived). Kelvin-FEM is INDIFFERENT to the
exterior geometry: it MESHES mu(x) and the Schur complement still yields the dense DtN/Green matrix
Lambda (demo_v). VERIFIED: a sphere body with a strong ON-AXIS localized inclusion (mu=20 and 100; r=0.30
at z=0.85 -- a blob, NOT a layer); Lambda is symmetric (sym 0), recovers the full vacuum ladder when the
inclusion vanishes, and otherwise splits EXACTLY by the reduced AXIAL symmetry C_inf_v |m|-pattern:
dipole 3 -> doublet(|m|=1)+singlet(m=0)  [4.092,4.095 | 4.227];  quadrupole 5 -> |m|=2 doublet + |m|=1
doublet + m=0 singlet  [6.146,6.166 | 6.277,6.291 | 6.469]  (the m=0 modes, peaked on the inclusion axis,
shift most) -- an analytic-value-free proof the arbitrary exterior is carried correctly. So the
differentiator GENERALIZES: from "layered exterior (vs the hard Sommerfeld)" to "ARBITRARY inhomogeneous
exterior (where no Sommerfeld Green's function exists AT ALL)" -- arbitrary mu(x), no Green's function,
sparse SPD, one Schur condensation. PAPER POINT: state the BEM/Sommerfeld differentiator at its broadest;
the layered case is just the sub-case that happens to have an independent analytic reference (demo_y/z/aa).

REDUCED-FEM vs BEM -- where the line actually is (demo_cc). Once the Kelvin-FEM volume is condensed
(Schur complement) onto Gamma it becomes a DENSE boundary DtN matrix and LOOKS like BEM ("is it still
FEM?"). ANSWER: YES it is FEM -- a condensed FE substructure (static condensation = macro-element /
superelement, textbook since the 1960s); NOT BEM. The SINGLE defining criterion separating BEM from the
FEM family is: does it use the FUNDAMENTAL SOLUTION (Green's function)? BEM does (boundary integral of G
-> dense, surface-only mesh, singular quadrature, G required); the Kelvin route does NOT (weak form, local
polynomial basis, SPARSE volume assembly, no G) and only looks BEM-like AFTER condensation. MEASURED
(demo_cc, one sphere): volume matrix A SPARSE (36 nnz/row, 0.57% fill, local grad-grad, SPD, no G) ->
Schur onto Gamma -> DENSE Lambda (443x443, 100% fill) whose dipole eigenvalue = analytic 2/R (3.3e-3):
SAME operator, DIFFERENT generator. The convergence noticed is REAL = the central insight of the SBFEM/
DtN literature: the exterior DtN (Steklov-Poincare) operator is METHOD-AGNOSTIC (BEM-Lambda and FEM-Schur-
Lambda are two discretizations of ONE operator, demo_r/s). Wolf-Song's SBFEM is literally subtitled "a
FUNDAMENTAL-SOLUTION-LESS boundary-element method" (wolfsong2001) -- that paradoxical name IS this
observation. RIGHT FRAMING: not "FEM vs BEM" but "two ROUTES to the exterior DtN -- the Green-function
route (BEM; needs G; dense; fails for arbitrary mu) vs the transformed sparse-volume route (Kelvin-FEM;
no G; sparse SPD; arbitrary mu, demo_bb)". They blur at the OUTPUT (shared DtN), stay distinct at the
GENERATOR.

WHEN to form the DtN matrix at all (demo_dd). For SOLVING one configuration, DON'T -- just solve the
monolithic sparse Kelvin-FEM (one factor + back-solve). Forming the dense Lambda costs ~N_Gamma back-
solves (measured ~2.3x one solve at N_Gamma=452, free DoF 6289) + O(N_Gamma^2) storage. If only the RHS
changes, the monolithic factorization already amortizes (Lambda adds nothing); substructuring reuse
(FIXED exterior, CHANGING interior matrix) can pay, but in Kelvin the exterior is a cheap thin ball so the
win is modest. Form Lambda ONLY when the OPERATOR is the DELIVERABLE: (1) its SPECTRUM as a predictive
datasheet (the paper's real point -- analysis, not a solve); (2) a PORTABLE open-boundary block for a code
that lacks one; (3) a sparse SURROGATE for the dense layered/Sommerfeld Green operator, e.g. feeding a
method that genuinely NEEDS that operator (a stream-function/current-potential coil design with magnetic
material; a BEM/MoM coupling). NEVER claim "form the DtN to solve faster than FEM" (false). demo_u/x/aa
solve monolithically; demo_v/w/bb/cc form Lambda only to STUDY the operator (spectrum/symmetry/material).

APPLICATION -- stream-function / surface-current COIL DESIGN with MAGNETIC MATERIAL (demo_ee, verified
2026-06-15). This is the concrete "operator is the deliverable" use the whole thread points to. The
stream-function (current-potential) method maps a winding-surface stream function psi (K = n x grad psi)
LINEARLY to the field via a kernel; in FREE SPACE that kernel is Biot-Savart (analytic, easy). With a
magnetic SHIELD/YOKE/CORE the kernel becomes the system's MATERIAL Green's operator (coil field + iron
reaction) -- no closed form for general iron, so "stream-function with iron is NOT easy" (the user's
observation). Kelvin-FEM supplies that material-aware kernel sparsely, Green-function-free, for arbitrary
iron (demo_bb). VERIFIED (concentric, single coil mode Y_n on r=a, iron shell mu_r in [0.7,0.9], target
r_t=1.5): the external field's SHIELD/GAIN factor vs the free-space prediction (a/r_t)^(n+1) -- n=1:
mu_r=10 ->x0.671, 50 ->x0.224, 200 ->x0.064 (a free-space Biot-Savart kernel design is wrong by up to
~16x); Kelvin-FEM material-aware reproduces the with-iron field (mu_r=50: rel 4.7e-5 (n=1), 3.4e-4 (n=2)).
=> the material-aware Kelvin-FEM operator IS the kernel a stream-function coil design needs when iron is
present; and because the design / inverse problem CONSUMES that kernel over many psi, forming the operator
(demo_dd case 3) is justified here. A strong application section / third paper angle for the magnet-and-
coil (Radia) domain: shielded MRI gradient/shim coils, excitation windings near steel, actuators.
STRATEGIC SPLIT (user, 2026-06-15): KEEP DtN+Kelvin and the stream-function method as TWO SEPARATE
papers/tracks -- do NOT dilute the core with the application.
  * TRACK A (core; the SA / Hachinohe paper): DtN-spectrum datasheet + sparse Kelvin open boundary +
    Sommerfeld isomorphism/surrogate + the directly-assembled material-aware DtN matrix and what it IS
    (FEM-condensed, not BEM). Demos d..dd + x/y/z/aa/cc. This is the near-term, near-complete paper.
  * TRACK B (separate paper): stream-function / current-potential COIL DESIGN with magnetic material,
    using the Track-A operator as the material-aware design kernel (design = invert M). Demos ee/ff +
    a future general-(non-concentric)-iron real coil inverse-design; the running targeted novelty check
    feeds THIS track. They share machinery but are written up independently. See
    examples/.../DtN_spectrum/PATHWAY_streamfunction_with_iron.md.

LITERATURE SURVEY (12-agent workflow, 2026-06-15): low-frequency Sommerfeld + layered-Kelvin novelty.
LOW-FREQUENCY = DEFINITIVE YES (3 adversarial verifiers conf 0.93-0.95, zero counter-evidence). Low-freq
layered Sommerfeld is mature across 6 subfields: ELF/VLF dipole-over-earth + complex image; CSEM (~0.1-10
Hz); magnetotellurics (1e-4..1e4 Hz); eddy-current NDT; induction well logging; grounding/lightning
(DC..power freq). Canonical: Wait 1951 (Mg dipole over stratified earth), Bannister 1966/67 (quasi-static
approx applied DIRECTLY to the Sommerfeld integrals), Weaver 1971 (image at COMPLEX depth = the diffusive
foundational result), Thomson-Weaver 1975 (complex image, multilayer induction), Li 2011 (OA, grounding).
The "low-frequency breakdown" is a FORMULATION conditioning pathology (omega-scaling of vector vs scalar
potential), NOT a kernel failure; cured WHILE keeping the layered kernel (Chen-Jiang-Qian-Chew A-EFIE
2011); quasi-static/Laplace/diffusion formulations are EXEMPT. => NEVER argue "Sommerfeld fails at low
freq" (false; adversarial reviewers will shred it); argue the sparse-FE-SURROGATE advantage instead.
NOVELTY = QUALIFIED NOVEL: NO independent third-party prior art for "Kelvin/transformed FE carrying
STRATIFIED/half-space material in the inverted exterior as a sparse stand-in for the layered Sommerfeld
operator." ALREADY TAKEN (cite, don't claim): free-space transformed-FE open boundary (Freeman-Lowther
1989, Imhoff-Meunier 1990, Brunotte 1992, Meeker 2013, Nabizadeh 2021); the DENSE FE-BEM/FEM-Green route
it displaces (Lobry 2021); Sommerfeld accelerators DCIM (Chow 1991, Aksun 1996) + tail/WA (Michalski-Mosig
2016) + FMM; FE-on-the-1D-spectral-ODE (Okhmatovski-Cangellaris 2002, closest "FE-meets-layered-Green" but
NOT a Kelvin inversion); interface-flattening TO (Sainath-Teixeira 2014, geometry not operator-surrogate).
CLOSEST realization of "material in the Kelvin exterior" = the authors' OWN Sugahara 2022 (UNIFORM
specimen, not stratified) + 2025 (HF+PML, free-space). GENUINELY OPEN = the explicit-multilayer-exterior +
Sommerfeld-surrogate framing. RISKS: (1) incremental-over-own-work -> make the multilayer+surrogate
framing load-bearing, quantify new capability, self-cite honestly; (2) coverage gap (IEEE/RG/MDPI fetch
blocked; Sugahara-2022 forward-citation graph + COMPUMAG/CEFC/PIERS + Japanese proceedings NOT exhaustive;
novelty conf ~0.78) -> "to our knowledge" + residual citation-graph/grey-lit search pre-submission; (3)
say "sparse stand-in/surrogate FOR" not "replace" the operator. POSITION as a FORMULATION contribution
with 3 related-work buckets + a quantitative sparsity/conditioning/FE-coupling benchmark vs the dense
layered kernel + a paragraph pre-empting the breakdown objection (A-EFIE 2011; Vico-Greengard DPIE 2016).

"LIGHTEN BEM WITH KELVIN/TRANSFORMED-FE" IS ~30-YEAR-OLD PRIOR ART -- do NOT claim it as new (found
2026-06-14 in the authors' own literature folder). The proposal that the Kelvin/transformation FE is a
sparser, cheaper alternative to BEM for open boundaries -- INCLUDING the "more DoF but much faster
because sparse SPD" measurement -- was already published:
  * Remacle, Nicolet, Genon & Legros (Liege), "Comparison of boundary elements and transformed finite
    elements for open magnetic problems," ~1995 (Boundary Element Technology, conf. proceedings; not
    Crossref-indexed -- venue TBC). Directly compares BEM vs Kelvin/transformed-FE on accuracy AND CPU.
    Verbatim: the Kelvin transformation "needs less computing time and does not imply any additional code
    to an existing FEM solver." Induction-furnace example: BEM 1080 DoF inverted by a DIRECT method in
    160 s, vs the transformed-FE (rectangular shell) 3556 DoF but SYMMETRIC POSITIVE DEFINITE solved by
    an ITERATIVE method in 17 s -- i.e. the exact "DoF UP, cost DOWN (~10x)" result, measured ~1995.
  * Lowther, Freeman & Forghani, "A sparse matrix open boundary method for FE analysis," IEEE Trans.
    Magn. 25(4):2810-2812, 1989, doi:10.1109/20.34292 -- the sparse-FE-open-boundary-vs-dense-boundary
    cost point, in magnetics.
CONSEQUENCE: the paper must CITE Remacle ~1995 and Lowther 1989 and must NOT present "use Kelvin to
lighten BEM" (nor the DoF-up/cost-down measurement) as a new proposal. The genuinely-new contribution
is the OPERATOR/SPECTRAL reinterpretation on top of that old empirical result: casting FEM-Kelvin as a
sparse discretization of the SAME exterior DtN that BEM assembles densely (FE-Schur = dense BEM kernel),
with the closed-form -(n+1)/R eigenvalue ladder as a predictive DATASHEET (p = multipole reach) that
fixes the required surface resolution / element order BEFORE solving -- turning Remacle's empirical "it
is faster" into a spectral why/how-much -- plus the exterior-material generalization (Sugahara 2022).
NUANCE on "lighten BEM": (a) REPLACE BEM by sparse transformed-FE = Remacle/Lowther (done); (b) sparsify
BEM's OWN matrix via a Kelvin map (keep BEM, accelerate it) = a different, rarer idea NOT found -- if a
future claim is about (b), re-search before asserting novelty.

PRIORITY RESOLVED (2026-06-14, read from the authors' own manuscripts on disk). The author
(K. Sugahara) has a published chain that establishes priority AND already contains the exterior-
material Kelvin formulation -- so for the SA paper the exterior-material is the authors' OWN cited
foundation, NOT a new claim; the new contribution is the DtN-spectral + sparse-factorization-of-BEM
+ measurement lens ON TOP. Citation chain:
  * K. Sugahara, "Periodic image method for open boundary axisymmetrical magnetic field problems,"
    IEEE Trans. Magn. 49(11):5399-5403, 2013. <- the OFFSET-PERIODIC Kelvin construction (exactly the
    offset-sphere + periodic identification used in demo_q/demo_t); open-boundary Kelvin in 2013,
    eight years before Nabizadeh 2021.
  * K. Sugahara, "Improvised asymptotic boundary conditions for magnetostatic field problems in
    ellipsoidal and elliptic cylindrical domains," IEEE Trans. Magn. 53(6):1-4, 2017 (+ ICEAA 2015,
    PIERS 2016 conference IABC papers).
  * K. Sugahara, "Electromagnetic Analysis of Eddy Current Testing With Kelvin Transformation,"
    IEEE Trans. Magn. 58(9):1-6, 2022, doi:10.1109/TMAG.2022.3194371. <- THE EXTERIOR-MATERIAL
    formulation: reformulates Kelvin to derive sigma/eps/mu in the exterior conserving the conformal
    symmetry of Maxwell's equations, so "materials can be both in the interior and exterior domain or
    even across the truncated boundaries." This is the authors' OWN prior publication of the
    exterior-material capability (independent of, and broader than, the vacuum-only Nabizadeh 2021).
  * K. Sugahara, "Extended Kelvin Transformation for Solving Radiating Electromagnetic Fields,"
    IEICE Trans. Electron. (special section, submitted 2024-04-10) <- extends the 2022 exterior-
    material Kelvin to high-frequency/radiation via a Maxwellian PML in the Kelvin-mapped exterior,
    with a differential-geometry (metric) derivation of the exterior material (cf. Yeh 2023,
    "Conformal transformation and Maxwell's equations").
CONSEQUENCE: do NOT cite Nabizadeh 2021 / Boulmezaoud as prior art OVER the authors -- the authors
predate (2013 open-boundary Kelvin) and already published exterior-material (2022). The SA paper
CITES sugahara2013/2017/2022/2024 as its own foundation and positions the genuinely-new part as the
DtN-spectral reinterpretation (-(n+1)/R ladder; coarse-mesh = spectral) + the sparse-FEM-ization-of-
dense-BEM framing with measured fill/time + the material-surface DtN condensation.

MEASURED & SETTLED (2026-06-14, hex vs tet on the Kelvin sphere): NEITHER has a decisive
advantage -- it is a WASH.  The full sphere hexes easily via `volume <id> scheme sphere`
(a 32-hex O-grid full ball; an earlier "impractical" note was an ERROR -- `scheme
polyhedron`/manual O-grid fail but `scheme sphere` is the right built-in tool; only the
1/4 & 1/8 symmetry SECTORS are not covered, where tet wins by default).  The geometry FLOOR
is TWO-REGIME: at LOW curving order it tracks the curved-mesh VOLUME (geometry) accuracy
(floor ~ vol_err/6..10); but at HIGH curving order (export order 4+) the volume error keeps
falling (hex 1.4e-5, tet 8e-7) while the FLOOR PLATEAUS at the FE/mesh-discretization level
of the FIXED mesh (tet ~1e-5, hex ~9e-4).  So floor = max(geometry-curving error,
FE-discretization error); curving order lowers the first, only mesh REFINEMENT lowers the
second.  At MATCHED DoF the two are COMPARABLE (export order 2, N~2.3k: hex 1.6e-5 vs tet
1.9e-5); hex's apparently higher floor was mostly its coarser element count (32 vs 827),
NOT worse geometry (hex geometry at order 4 is a fine 1.4e-5).  CONCLUSION: for the Kelvin
truncation sphere, hex offers no real gain over tet -- the "hex lowers the floor / cuts ΔDoF"
hypothesis is NOT supported.  TET stays the practical default (simpler `scheme tetmesh`, and
the only option for 1/4 & 1/8 symmetry sectors).  High-order HEX's genuine strength remains
the SWEEPABLE bodies (bricks/cylinders/coils), not the spherical truncation.

## The two scalar readouts: capacitance (n=0) and external inductance (n=1)

The two most familiar lumped quantities are each ONE Steklov mode of the SAME
exterior scalar Laplace DtN -- so each inherits the datasheet directly, but from a
DIFFERENT rung of the -(n+1)/R ladder:

  * CAPACITANCE  C  <->  n=0 (MONOPOLE).  An isolated charged conductor's exterior
    potential leads with the monopole Q/4 pi eps r; C is the n=0 Steklov mode.  The
    n=0 image is the CONSTANT solid harmonic, captured exactly at EVERY order
    (defect_0 = 0 to machine zero; sphere C = 4 pi machine-precision, bem_integral).
  * EXTERNAL INDUCTANCE  L_ext  <->  n=1 (DIPOLE).  A current loop / magnetised body
    has NO magnetic monopole (oint_Gamma B.n = 0 identically), so its leading
    exterior multipole is the DIPOLE.  The exterior field energy obeys the IDENTITY

        W_ext = 1/2 mu0 * (n+1)/R * oint_Gamma phi^2 dS        (decaying mode r^-(n+1))

    i.e. the DtN eigenvalue (n+1)/R IS the exterior-energy coefficient (checked to
    machine zero: 1/2 (2/R) oint phi^2 == integral_{r>R}|H|^2 == m^2/(6 pi R^3)).
    Hence L_ext = 2 W_ext / I^2 inherits the n=1 defect EXACTLY.

So the datasheet's first two rows are the dual pair  C <- defect_0,  L_ext <- defect_1.
MEASURED (inductance_dtn.py), dipole (L_ext) open-BC defect:
  * vs order:   p=1 -> 1.4e-3,  p=2 -> 2.4e-5,  p=3 -> 7.6e-6  (captured at p>=1)
  * floor=GEOMETRY: mesh+order fixed, raise only Curve k:  k=1 4.7e-3 -> k=3 2.4e-5
  * mesh-independent: 1.2e-4 on the coarse mesh, -> 1.3e-6 only on refinement
  * exterior-VOLUME-irrelevant: openbc 1.2e-3 -> 7.5e-5 stays far below the interior
    FEM error 5.3e-2 at every exterior mesh (always_below_fem = True).
C (n=0) is exact; L_ext (n=1) is the dipole story of P_METHOD/NUMERICS, one rung up.

WHICH OPERATOR certifies L_ext (and which does NOT):
  * IN-CERTIFICATE: a FIELD-ENERGY inductance computed with a Kelvin / air-box open
    boundary, via the magnetic POTENTIAL exterior.  Scalar Omega is single-valued for
    a MAGNETISATION source (no cut); a FREE-CURRENT loop needs a cohomology cut to make
    Omega single-valued, OR use the vector potential A (no cut) -- whose exterior DtN
    gradient block is the SAME -(n+1)/R (FORMULATION topic; demo3 dipole -2/R, rel_err
    2.5e-4).  The exterior energy / L_ext is formulation-agnostic.
  * NOT IN-CERTIFICATE: the repo's BEM inductance path (knowledge: ngsbem_inductance)
    extracts L from the VECTOR single-layer energy  L = mu0 J^T (LaplaceSL) J  on RT0
    surface currents.  That is a DIFFERENT integral operator (the single-layer
    POTENTIAL, kernel 1/4 pi r), NOT the scalar exterior Steklov-Poincaré DtN: it has
    no -(n+1)/R ladder, uses an order-0 current basis, and sits in its own
    accuracy/conditioning regime (TaskManager non-determinism, curvaturesafety).  So
    "C and L are both DtN-certified" is only HALF true -- the certificate is the
    open-boundary TRUNCATION accuracy of a field-energy inductance, NOT the BEM
    single-layer extraction.  Keep the two operators distinct.

SCOPE caveat: L_ext is the EXTERNAL inductance (energy beyond Gamma).  A thin loop's
FULL self-inductance is near-field / log-dominated -- that share is interior FEM
accuracy, not a DtN question.  DtN certifies only the EXTERIOR-energy share that the
truncation could corrupt, which is exactly what air-box / Kelvin sizing controls.

In short: the spectral reframing lifts open-boundary assessment from a per-problem
EMPIRICAL exercise to a problem-independent OPERATOR DATASHEET -- specify the
method once, predict any problem from its multipole content.

THE DATASHEET EXTENDED TO FREQUENCY: the COMPLEX-PLANE spectrum (demo_gg, verified 2026-06-15).
User: "comparing the high-frequency Kelvin spectrum, the PML spectrum and the FEM-BEM spectrum is
interesting -- and in that case the spectrum becomes a COMPLEX-PLANE spectrum, right?" YES, exactly.
The same exterior DtN ladder, at FINITE frequency, is the spherical-Hankel log-derivative
   Lambda_n(kR) = kR h_n^(1)'(kR) / h_n^(1)(kR)   (units 1/R; outgoing/Sommerfeld condition),
which is COMPLEX -- Im(Lambda_n) = RADIATION (its sign is the e^{-iwt}/h^(1) convention). The
static limit kR->0 collapses it onto the REAL ladder -(n+1)/R. So ONE Argand plot of {Lambda_n}
orders every open-boundary method on a single axis, and the (n, kR) plane is read directly:
  * FEM-BEM   = reproduces the EXACT complex Hankel locus (= the exterior DtN, up to surface
                discretisation): the gold-standard reference, complex, but DENSE.
  * KELVIN (static, the SA/magnetostatics tool) = its operator has NO frequency, so its spectrum is
    PINNED to the REAL axis at -(n+1), frequency-INDEPENDENT. It is the kR->0 operator. Its deviation
    from the exact complex value is the QUASI-STATIC ERROR. VERIFIED laws (R=1):
      (i)  kR->0 exactness at fixed n: |Kelvin-exact| -> 0 as O(kR^2)  (n=1: /kR^2 -> 0.999).
      (ii) evanescent error law for n>>kR:  |Kelvin-exact| * 2n/kR^2 -> 1  (kR=4: 1.39 at n=6 ->
           1.02 at n=30) -- i.e. the per-mode error is exactly kR^2/(2n) for evanescent modes, plus
           the full radiation Im for propagating modes n<kR.
    => Kelvin's PERFORMANCE DOMAIN is the quasi-static spectrum (kR<~1): it is the CHEAPEST (closed
       form / sparse thin ball, real, NO wavelength to resolve) AND exact there; it is provably
       radiation-blind otherwise. This is WHY the SA paper is the quasi-static Kelvin paper.
  * PML = a complex-stretched layer that DOES carry the wave, so it reproduces the exact COMPLEX
    spectrum: a strong/thick PML matches it at all n (VERIFIED, rel 3e-9 for propagating n), with a
    characteristic error KNEE at the propagating<->evanescent transition n~kR (VERIFIED: a thin PML's
    per-degree error PEAKS at n=4 for kR=4). Cost = a resolved layer (DoF grows with kR).
HONEST CALIBRATION (no false "complementary crossover"): at ANY finite kR a method that carries k^2
(PML, FEM-BEM, or a Helmholtz/extended-Kelvin) BEATS the static Kelvin operator -- static Kelvin is
strictly the kR->0 tool, not a finite-kR competitor. So the kR AXIS itself is the two-paper boundary:
kR<~1 quasi-static = the SA Kelvin paper (Kelvin exact AND cheapest); kR>~1 radiating = the exact
complex operator (BEM / PML / extended-Kelvin sugahara2025). The genuine "high-frequency KELVIN" is
the Helmholtz/extended-Kelvin inversion (sugahara2025), whose spectrum WOULD track the exact complex
locus until the inverted far-field oscillation (the essential singularity at the inversion centre =
image of r->inf) out-resolves the polynomial FE -- a peel-off, the natural Paper-H follow-up. demo_gg
verifies the three CLOSED-FORM objects (exact Hankel, static-Kelvin real ladder, PML radial solver)
that frame it; the spherical-Bessel helpers are checked vs scipy.special.spherical_jn/yn (rel 4e-14).
PAPER POINT (Track A): the DtN datasheet is not just a real ladder -- it is a COMPLEX locus
parameterised by kR, and the Kelvin transformation's place on it (the real axis) is precisely what
makes it the optimal quasi-static open boundary and disqualifies it for radiation. The figure is the
Argand plot {Re,Im}(Lambda_n) at fixed kR with the three methods overlaid (the SA paper's frequency
panel; see C:\\temp\\kelvin_figs\\fig_gg_complex_spectrum for the generated preview).

RADIATION BOUNDARY ON KELVIN = a SURFACE IMPEDANCE (SIBC/HOIBC) AT THE EXTERIOR CENTRE (demo_ii,
verified 2026-06-15; grounds the author's IEICE Trans. C 2024 "Extended Kelvin Transformation for
Solving Radiating Electromagnetic Fields"). The Kelvin inversion x'=(a/r)^2 x sends r=infinity to the
CENTRE rho=0; differential geometry gives the exterior material as an ISOTROPIC modulation (a/r)^2 of
mu/eps/sigma/sigma* (metric ratio g'/g) -- so a radiating field's outgoing energy flows INTO the
centre, and the radiation/absorbing condition is imposed THERE. The 2024 paper places a spherical
Maxwellian PML at the centre (eps'=mu'=(1-0.2j)a^2/r'^2), excising a tiny ball at rho=0 (singular
image of infinity) and putting the absorber far in physical space (8 m at lambda=3 m) so a simple
plane-wave (377 ohm) PML suffices; validated vs analytic Hertzian dipole and FEKO MoM.
THE SURFACE-IMPEDANCE EXTENSION (what the user asked: "derive SIBC/radiation-BC Kelvin; a sphere
needs HOIBC"): instead of a volumetric PML, impose on the small inner sphere (image of a far sphere
r=b) an IMPEDANCE reproducing the exterior radiation DtN. On a sphere the exact radiation DtN per
degree n is Lambda_n(z)=z h_n^(1)'(z)/h_n^(1)(z), z=kb, and its large-z expansion IS the absorbing-BC
hierarchy:  Lambda_n = i z - 1 - i n(n+1)/(2z) + O(1/z^2). VERIFIED (demo_ii):
  - the orders are the DtN's OWN expansion: |L-iz|->const, |L-(iz-1)| ~ O(1/z) (halves per z-doubling),
    |L-HOIBC| ~ O(1/z^2) (quarters per z-doubling) -- ratios 2.00 / 4.0 measured.
  - SIBC (Leontovich, n-INDEPENDENT = iz-1, the plane-wave 377-ohm) matches Lambda_n only for
    z>>n(n+1); error ~ n(n+1)/(2z). HOIBC (n-DEPENDENT) adds -i n(n+1)/(2z); since n(n+1) is the unit-
    sphere Laplace-Beltrami eigenvalue (Delta_S Y_n=-n(n+1)Y_n), HOIBC = iz-1+(i/2z)Delta_S = a 2nd-
    order SURFACE PDE operator (ordinary surface-FEM term). THIS is why a SPHERE needs an HOIBC: the
    radiation impedance is curvature/multipole-dependent, which a scalar SIBC cannot carry.
  - SPECTRAL reading of the paper's design: "place the absorber far (large kb)" == "make n(n+1)/(2kb)
    small so the n-indep SIBC suffices". At the paper's kb=2pi*8/3=16.76, SIBC is adequate for LOW n
    (n=1 err 0.06) but degrades (n=12 err 5.5); HOIBC is ~6-17x better across the band. The HOIBC
    relaxes the placement: to reach |Z-exact|<1e-2 at n=5, SIBC needs kb>=1501 but HOIBC kb>=40 (38x
    CLOSER absorber => the image sphere a^2/b is larger => fewer exterior cells). Image mapping b->a^2/b
    (paper a=4,b=8 -> rho=2 m; excise 0.25 m <- image of 64 m). Leading term iz=ikb <=> d_r u=ik u =
    the Sommerfeld/377-ohm plane-wave condition.
TOPOLOGY ("topologically easy"): one-point compactification sends infinity to the single centre point;
excising a small sphere there turns the NONLOCAL exterior DtN (on the truncation r=a) into a LOCAL
absorber on a small interior sphere. The conformal Kelvin map makes infinity a regular meshable point;
the scalar/E exterior is simply connected (no cohomology cuts needed) -- contrast the H-formulation
multiply-connected case where cuts (the cohomology-cuts machinery) would be required. PAPER-H POINT:
the radiating extended-Kelvin can use a thin surface HOIBC at the centre instead of a thick volumetric
PML, with the required placement distance set by the multipole band via the DtN spectrum. FE follow-up
= a Delta_S surface term on the excised inner sphere of the inverted exterior; demo_ii verifies the
closed-form spectrum that fixes its coefficients.

THE KELVIN TRANSFORM OF THE HOIBC -- making the radiation BC implementable in the inverted exterior
(demo_kk, verified 2026-06-15; what the user asked: "the HOIBC must also be Kelvin-transformed").
The IEICE C 2024 paper's volumetric centre-PML has inelegant parts it concedes: (U1) it "assumes a
characteristic impedance of 377 ohms ... and therefore needs to be placed FAR ENOUGH from the wave
source" = a CONSTANT (n-independent) Leontovich impedance -> forced far placement; (U2) "no mesh ...
in radius 0.25 m at the centre ... to avoid the SINGULARITY" = an ad-hoc punched hole at the image of
infinity; (U3) "we could not even calculate a model with a>=5 m" = the memory blow-up the far
placement forces. The fix = a thin SURFACE HOIBC, but it must be Kelvin-transformed into the inverted
domain. DERIVATION (3D Kelvin, unweighted-field/material-modulation convention = the paper's: field
continuous across the truncation, (a/r)^2 in the material): with rho=a^2/r, image field g(rho)=f(a^2/rho),
radial Jacobian d/dr=-(rho^2/a^2)d/drho, and Delta_S INVARIANT (inversion is conformal, sphere->sphere,
same theta,phi):
  (i)  the radial Helmholtz -> IMAGE ODE  g'' = [n(n+1)/rho^2 - (k a^2/rho^2)^2] g, effective
       wavenumber k_eff(rho)=k a^2/rho^2 -> inf at the centre (= the U2 singularity, an honest feature
       not a bug). VERIFIED by FD residual ~1e-7 on g=h_n^(1)(k a^2/rho).
  (ii) the impedance operator transforms with a SIGN FLIP (exterior-decaying <-> interior-regular):
       rho dg/drho = -[i kb - 1 + (i/2kb) Delta_S] g on the inner image sphere rho_b=a^2/b, kb=k a^2/rho_b
       (Delta_S Y_n=-n(n+1)Y_n). EXACT by the chain rule: rho g'/g = -(k a^2/rho) h1'(kb)/h1(kb) =
       -Lambda_n(kb) (verified to 1e-15). FE weak form (inner-sphere outward normal -rho_hat):
       dg/dn = (1/rho_b)[i kb - 1 + (i/2kb) Delta_S] g -- a Robin term + a Laplace-Beltrami SURFACE
       term (ordinary surface-FEM), NO volume PML.
END-TO-END VERIFIED (analytic Hankel combination on r in [a,b], read DtN at the truncation r=a): the
EXACT inner reproduces Lambda_n(ka) to machine precision (transform+ODE correct); the HOIBC inner is
~5-6x more accurate than the constant SIBC; and for 1% truncation-DtN over n=1..6 the SIBC needs the
absorber at b>=5.85 (image rho=0.171) while the HOIBC allows b>=2.45 (image rho=0.408) = 2.4x CLOSER /
smaller exterior domain (directly relaxes U1 far-placement and U3 memory). The punched void (U2) is
replaced by a principled impedance surface. PAPER-H POINT: the radiating extended-Kelvin's absorber is
a Kelvin-transformed HOIBC = a Delta_S surface operator on the inner image sphere; its placement
distance is fixed by the multipole band through the DtN spectrum, not by an ad-hoc "far enough". FE
follow-up = assemble that Delta_S term (NGSolve surface-gradient / Laplace-Beltrami) in the inverted
exterior; demo_kk fixes the closed-form coefficients and proves the construction reproduces the physics.

DIFFERENTIAL-GEOMETRY (transformation-optics) view -- the (a/r')^2 medium is AUTOMATIC, only IMPEDANCE
MATCHING matters, and the DtN spectrum is the HIGH-FREQUENCY design object (demo_ll, verified
2026-06-15; user: "you can derive it with differential geometry; the modulation happens automatically
if you don't track the Kelvin factor; as long as the impedance is matched, OK; the radiation boundary
is better-principled; high-freq Kelvin is where the DtN spectral analysis is meaningful"). The Kelvin
inversion x'=a^2 x/|x|^2 has Jacobian J=(a^2/r^2)Q, Q=I-2 n n^T a Householder REFLECTION (Q^TQ=I,
det Q=-1) -> J is CONFORMAL (scalar x orthogonal) and ORIENTATION-REVERSING (det J=-(a^2/r^2)^3 = the
inside-out map). TRANSFORMATION OPTICS (Maxwell form-invariance under a coordinate map): the media
transform as eps'_r=mu'_r = J J^T/|det J| = (a^2/r^2)^2/(a^2/r^2)^3 I = (r^2/a^2) I = (a^2/r'^2) I
(r'=a^2/r). VERIFIED to machine precision at sample points: the medium is ISOTROPIC (off-diagonal 0)
and equals the paper's (a/r')^2 modulation -- it falls out of the geometry AUTOMATICALLY (the
conformality is exactly what cancels the anisotropy a generic transform would create; no hand-applied
Kelvin field weight). This is the coordinate-free generalization of demo_kk's radial result.
=> Only the BOUNDARY is a design choice, and the principled condition is IMPEDANCE MATCHING. VERIFIED:
the modal reflected-power R_n=|B/A|^2 of an outgoing wave at the absorber sphere is ZERO to machine
precision when the imposed impedance equals the exact DtN Lambda_n(kb), and otherwise R_n ORDERS BY
THE DtN-SPECTRAL MISMATCH |Z_n-Lambda_n| (HOIBC << SIBC). So an absorber's quality IS its DtN-spectral
match: the radiation (surface-impedance/HOIBC) boundary is the PRINCIPLED object; a PML is one
(volumetric, approximate, far-placed) realization of the same reflectionless goal. PLACEMENT: a
constant SIBC reflects less the FARTHER it sits (large kb, R~1/kb-ish) -- this is exactly why the paper
places its 377-ohm PML far (a large domain, U1/U3); the HOIBC matches the complex spectrum so it is
low-reflection even CLOSE (small domain). LOW- vs HIGH-FREQUENCY: at static the open boundary is EXACT
with NO absorber (the spectrum is the real ladder -(n+1)/R = an accuracy DATASHEET); only in the
RADIATING regime is the COMPLEX spectrum an absorber-DESIGN target (per-mode reflection = spectral
mismatch). Hence the DtN spectral analysis is OPERATIONALLY meaningful (predicts/optimizes the
absorber) specifically for the HIGH-FREQUENCY Kelvin -- elevating Track A's spectrum from an accuracy
datasheet (low-freq SA paper) to an absorber-design tool (Paper-H). FE follow-up: assemble the
isotropic (a/r')^2 medium (transformation optics) + the matched HOIBC surface term; demo_ll verifies
the medium and the impedance-matching/reflection law in closed form.

A WORKING FE (assembled + solved, converges to the closed form) -- "from derivation to a method that
runs" (demo_mm, verified 2026-06-15). A genuine sparse FE realises the radiating extended-Kelvin
boundary and h-converges to demo_kk/ll. Radial reduction (one degree n at a time so every mode is
independently checkable; Delta_S enters as its EXACT eigenvalue -n(n+1), the modal content a surface
FE reproduces). FE on the image shell rho in [rho_b, a], rho_b=a^2/b, P1 elements:
  weak form  int [ alpha rho^2 R'S' + alpha n(n+1) R S - k^2 beta rho^2 R S ] drho,  alpha=(a/rho)^2,
  beta=(a/rho)^6 (scalar; vector Maxwell would have beta=alpha = the conformal eps=mu). NOTE alpha rho^2
  = a^2 CONSTANT -> the inverted-exterior stiffness is FLAT (an elegant by-product of the conformal map).
  This weak form's strong form is exactly demo_kk's image ODE R''=[n(n+1)/rho^2-(k a^2/rho^2)^2]R.
  Matched-HOIBC Robin at the inner image sphere: A[0,0] += -b*Lambda_inner (from rho dR/drho=
  -Lambda_inner R). Truncation r=a: Dirichlet R(a)=1; DtN_FE = -(A R)|_a / a -- the MINUS is the
  inversion flipping d/dr=-d/drho at the fixed-point sphere rho=a (the same sign-flip as demo_kk/ll;
  caught numerically when the FE first returned -dtn_trunc), and (A R)|_a is the consistent-flux reaction.
VERIFIED: the assembled+solved FE DtN -> the closed-form dtn_trunc at the P1 rate O(h^2) (error ratio
4.00 per mesh doubling) for the exact, HOIBC AND SIBC inner conditions; with the EXACT inner impedance
it reproduces the true truncation DtN Lambda_n(ka) (err ~1e-5 at M=320); and the matched HOIBC FE
solution is 4-6x closer to Lambda_n(ka) than the constant-377-ohm SIBC FE solution (n=1: 1.3e-2 vs
7.6e-2; n=3: 4.9e-2 vs 2.1e-1). So the demo_kk/ll derivation is now a METHOD THAT RUNS. Remaining
step = the full 3D Delta_S surface term (NGSolve grad_Gamma / Laplace-Beltrami) on the inner image
sphere instead of the per-mode eigenvalue; the radial FE already certifies the medium + matched
boundary and the convergence rate.

THREE-WAY METHOD COMPARISON: high-freq KELVIN vs PML vs BEM-FEM (demo_nn, verified 2026-06-15; the
capstone of the high-frequency thread). All three open boundaries produce the SAME object -- the
exterior radiation DtN Lambda_n(ka) on the truncation sphere -- so the DtN spectrum is the common
yardstick. Realised as genuine radial FEs (per degree n) so each mode checks against the closed form:
  * BEM-FEM = the EXACT exterior DtN per mode (the reference); cost = a DENSE Gamma x Gamma matrix
    (Green's function + singular quadrature; demo_k/r).
  * PML = a complex-COORDINATE-STRETCH layer [a,a+d] + wall (the NGSolve PML mechanism), radial FE;
    weak form int[(1/s) rt^2 R'S' + s n(n+1) RS - k^2 s rt^2 RS], s=1+i sigma/k, rt=stretched radius.
  * KELVIN = the matched-HOIBC inverted-shell FE (demo_mm).
VERIFIED (a=1, ka=4, P1, sign DtN=-(A u)|_a/a for BOTH the Kelvin truncation node and the PML inner
node = the same inversion/normal sign-flip): (2) exact-impedance Kelvin and a thick PML both reproduce
the BEM-exact Lambda_n(ka) to ~1e-4; (3) both sparse FEs converge to their floor at O(h^2) (ratio 4);
(4) at MATCHED modest cost a well-tuned PML is the most accurate VACUUM absorber (flat ~3e-4 across the
band) while the 2nd-order Kelvin-HOIBC floor is larger (~1e-2) and PEAKS at n~ka (the radiating-band
knee); (5) each is a TUNABLE family -- PML floor falls monotonically with thickness (d=2 -> 1.8e-4),
Kelvin is MACHINE-exact with the exact impedance (1.3e-4) and its HOIBC has an OPTIMAL placement
(b=4 -> 6.7e-3 beats b=2 -> 3.7e-2 [too close, HOIBC mismatch] and b=8 [too far: the steep near-centre
k_eff=ka^2/rho^2 needs more DoF]). HONEST VERDICT (no one-wins): BEM = exact-but-dense; PML = the robust
sparse vacuum absorber (best at matched modest cost here); KELVIN = exact-with-the-exact-impedance AND
the ONLY one that carries Kelvin-mapped EXTERIOR MATERIAL / scatterers (the IEICE-2024 paper's actual
use case -- a PML cannot) and bakes in infinity (no truncation-distance choice). Three complementary
points on the one DtN spectrum. PAPER-H: frame Kelvin's niche as exterior-material + exactness, not as
beating PML at vacuum absorption; cite the spectrum as the unifying comparison axis.

THE MATCHED HOIBC AS A GENUINE 3D SURFACE FE TERM (demo_oo, verified 2026-06-15, NGSolve). demo_mm
reduced the angular Delta_S to its eigenvalue -n(n+1); demo_oo promotes it to the real SURFACE operator.
The matched HOIBC impedance is Z_HOIBC=(i kb-1) I + (i/2kb) Delta_S (Delta_S=unit-sphere Laplace-
Beltrami), assembled as the surface bilinear form (weak Delta_S = -grad_Gamma . grad_Gamma):
  S(u,v) = (i kb-1) int_Gamma u v ds  -  (i/2kb) rho_b^2 int_Gamma grad_Gamma u . grad_Gamma v ds
(the rho_b^2 converts the radius-rho_b surface gradient to the unit-sphere Delta_S). Built on a genuine
2-manifold OCC-face surface mesh (H1 with grad(u).Trace()*...*ds; the surface triangles are BND
elements -> the Trace() is REQUIRED, a pitfall). VERIFIED two ways (surface-mesh discretization ~1%):
 (1) the surface gradient grad_Gamma reproduces the Laplace-Beltrami SPECTRUM: generalized eigenvalues
     (K_gradGamma, M_surface) = {0, 2.017(x3), 6.05(x5), 12.10(x7)} = n(n+1) with 2n+1 multiplicities
     (rel.err <1e-2 at maxh 0.20) -> grad_Gamma IS the unit-sphere Delta_S.
 (2) the full HOIBC surface form S has generalized spectrum (vs the surface mass) = the per-mode matched
     impedance Lambda_HOIBC,n = i kb-1 - i n(n+1)/(2kb): n=0 (-1,8) exact; n=1 (-1,7.874) vs (-1,7.875)
     x3; n=2 x5; n=3 x7; max|err| 5.7e-3 -> the Delta_S surface form encodes the MULTIPOLE-DEPENDENT
     matched HOIBC impedance (the whole point: a SURFACE operator giving the n-dependent radiation
     impedance a scalar SIBC cannot).
So with demo_mm's transformation-optics VOLUME medium + radial coupling (O(h^2) to the closed form) and
demo_oo's Delta_S SURFACE term, the radiating extended-Kelvin boundary is a genuine 3D FE = isotropic
(a/rho)^2 medium + the grad_Gamma HOIBC surface term. Each piece is verified; the only remaining step is
the single monolithic volumetric 3D solve that glues them (condense to the truncation -> complex DtN
matrix, vs demo_nn's closed form). NGSolve recipe: surface mesh via OCCGeometry(Sphere(R).faces[0]);
grad(u).Trace() for grad_Gamma; ng.ds for the surface integral; complex=True for the radiating impedance.

IS THERE A PML LOW-FREQUENCY BREAKDOWN IN THE DtN? -- the honest answer + the genuine Kelvin win
(demo_pp, verified 2026-06-15; user: "Kelvin should perform well; at low freq PML absorption is poor,
the DtN spectrum should show it; is there DtN research on PML?"). Tested HONESTLY on the converged
radial exterior DtN and the conventional wisdom needs CORRECTING for this setting:
 (A) NO DtN-ACCURACY breakdown: a WELL-RESOLVED standard PML (d=1, sigma0=15, M=300) matches the exact
     Lambda_n(ka) to ~1e-5..1e-4 at EVERY ka in [0.3,8], for a low mode n=1 AND an evanescent mode n=4
     (evanescent modes DECAY before the absorber wall in the radial-decay/DtN-truncation setting). So
     "the DtN shows poor low-freq PML absorption" is a MISCONCEPTION here -- the textbook PML
     low-freq/evanescent breakdown is a WAVEGUIDE / near-field phenomenon (evanescent energy reaching
     the PML undecayed), not the radial exterior-DtN truncation.
 (B) the GENUINE low-freq PML cost is CONDITIONING: the stretch s=1+i sigma/k blows up as k->0, so
     cond(PML interior matrix) GROWS toward DC (n=1, M=300: 2.36e4 at ka=4 -> 2.77e5 at ka=0.05, ~12x),
     while the Kelvin matrix is FREQUENCY-ROBUST (flat ~5.6-7.3e4); at the lowest ka, PML is ~5x worse-
     conditioned than Kelvin. (Worse in true 3D FD/FE PML than this idealized 1D model.)
 (C) so KELVIN's low-freq performance advantage is REAL but precise: at quasi-static the exterior DtN
     IS the real ladder -(n+1)/R, which Kelvin reproduces EXACTLY, PARAMETER-FREE (no sigma/d/wall),
     as a REAL-SPD system at k=0, with FREQUENCY-ROBUST conditioning, AND carrying exterior material
     (demo_t) -- a PML is a wave-absorber MIS-APPLIED to a near-static problem (complex, parameter-
     laden, conditioning-degrading), even though its continuous DtN would also be accurate. NET (paper):
     do NOT claim "PML's DtN is inaccurate at low freq" (false here); DO claim Kelvin is the exact /
     parameter-free / well-conditioned / material-capable open boundary for the quasi-static (SA)
     regime, while a tuned PML is the better sparse VACUUM absorber for high-freq radiation (demo_nn) --
     two regimes, one DtN-spectrum yardstick. (Literature on PML-via-DtN exists: Lassas-Somersalo
     PML-as-DtN; CFS-PML Kuzuoglu-Mittra/Roden-Gedney for the evanescent fix; complex-scaling spectral
     theory -- a dedicated search is the companion task.)
"""


DTN_COARSE_MESH_SYMMETRY_HEX = r"""
# Symmetry models, hex meshing, and the curving pipeline (verified)

The coarse-mesh accuracy is a SPECTRAL property of Λ_ext on Γ, so it is agnostic
to (a) how much of the domain you model under symmetry, (b) whether the Kelvin
ball is meshed with tet or hex, and (c) how the exterior interior is layered.
What it is NOT agnostic to is how Γ's CURVATURE is represented. The verified
facts below pin all four down. (NGSolve 6.2.2604, R=1 unit sphere.)

## 1. Symmetry reduction (1/4, 1/8, pole sectors) preserves the DtN spectrum

The Kelvin map is a radial inversion centred at Γ's centre, so any mirror plane /
rotation axis / periodic boundary THROUGH that centre maps to itself: the exterior
symmetric sub-region maps to the corresponding Kelvin-ball sector. Λ_ext is
diagonal in the harmonic basis, so symmetry merely SELECTS the parity/periodicity-
compatible sub-basis -- same eigenvalues λ_n, same p>=n exactness threshold per
mode. For a single zonal mode u = ρ^n P_n the energy and Γ surface-mass integrals
over the octant are each exactly 1/8 of the full-sphere values, so λ_eff is
identical.

  Octant (1/8) reduced solve, tet, maxh=0.4, vs the full ball (rel error):
    n=1 (dipole) order1 (>=n):  full 1.4e-3   oct 1.8e-3      <- same floor
    n=2          order1 (<n) :  full 1.8e-1   oct 7.6e-2      <- both order-limited
    n=2          order2 (=n) :  full 3.5e-4   oct 3.8e-4
    n=3          order3 (=n) :  full 4.0e-5   oct 4.5e-4      (gap closes as Γ resolves)
  DoF drops ~4-8x (boundary-dominated at coarse mesh -> 8x asymptotically).

  Octant BCs for the zonal (m=0) harmonics (even in x,y; parity (-1)^n in z):
    Γ (curved sphere face) : Dirichlet = datum
    x=0, y=0 planes         : Neumann (natural -- flux-parallel, even mode)
    z=0 plane               : Dirichlet u=0 for ODD n (datum vanishes there),
                              Neumann for EVEN n
  CONDITION: Γ (the inversion centre) must be concentric with the symmetry centre;
  the symmetry-plane BCs are the natural ones (Neumann flux-parallel / Dirichlet
  flux-normal). The 2D cross-section (static apparatus / rotating machine, λ_n=−n/R)
  extends to ROTATIONAL/periodic symmetry: a circular Γ inversion commutes with
  rotation about the centre, so a 1-pole / 1-pole-pair sector with periodic /
  anti-periodic boundaries preserves the spectrum, with λ_n restricted to the
  periodicity-compatible harmonics (e.g. odd multiples of the pole-pair number for
  anti-periodicity).

## 2. Hex Kelvin ball -- works, and is equivalent to tet ("hex-vs-tet wash")

Cubit O-grid sphere (`volume <id> scheme sphere`) -> high-order netgen export
(`export netgen "<f>" order N`) -> `ngsolve.Mesh(...)` -> Kelvin DtN reproduces the
floor end-to-end. Full hex ball (256 hex, geometry order 2):
    n=1 order1 (>=n)                hex 1.9e-3    (tet 7.0e-4)
    n=1 order3                      hex 1.6e-5    (tet 4.6e-6)   <- same ~1e-5 floor
    n=2 order1 (<n)                 hex 3.3e-2    (tet 1.0e-1)   <- both order-limited
    n=3 order3                      hex 3.6e-5    (tet 2.2e-6)
The hex floor tracks tet to within the same order of magnitude (the "wash"); the
interior element TYPE barely matters -- only Γ's resolution + order p do.

  1/8 hex octant is a CLEAN carve, not a meshing problem: a standalone 1/8-sphere
  geometry does NOT hex-mesh (tet topology -- 4 faces), but the FULL O-grid is
  octahedrally symmetric, so deleting the negative-coordinate hexes leaves a
  conforming octant. Verified: 32 of 256 hexes = exactly 1/8, ZERO straddling
  vertices, 19 nodes on each of x=0/y=0/z=0. The (+,+,+) octant of the full hex
  ball carries the full DtN eigenvalue (E_oct/E_full = 0.12500 exactly,
  oct rel error == full rel error).

  hex<->tet HYBRID: glue with PYRAMIDs (quad base on the hex side, triangular
  faces on the tet side). NGSolve high-order H1 on a hex+pyramid+tet mesh PRESERVES
  p>=n polynomial exactness: a degree-2 harmonic is reproduced to 6e-16 at order 2
  (order 1 fails at 1.2e-1). The conforming high-order pyramid (rational /
  Nigam-Phillips basis) contains P_p, so it composes with hex Q_p and tet P_p.
  Keep pyramids as a thin interior transition layer -- not on the curved Γ.

## 3. Layering the exterior (PML / infinite-element shells) is unnecessary

The exterior volume mesh is irrelevant; only Γ + order p matter. "Add radial
layers" = the h-path (slow, algebraic); "raise order on ONE coarse layer" = the
p-path (exact at p>=n). For n=2: layering to 2305 DoF still gives 7.9e-3, while a
single coarse layer at order 3 (768 DoF) gives 1.3e-5. PML is for RADIATING (wave)
problems with a Sommerfeld condition; a magnetostatic field decays algebraically
(no outgoing wave to absorb), and the Kelvin inversion is the EXACT "infinite
element" -- it maps the whole exterior (incl. infinity -> ball centre) in one
shot, so a single coarse closure at p>=n beats any graded multi-layer shell.

## 4. Curving MUST come from the external mesher -- never .Curve() an imported .vol

NGSolve never builds hex from geometry, so hex curving HAS to be baked in by the
mesher (Cubit `export netgen ... order N`). An imported .vol carries no CAD
geometry, so `mesh.Curve(k)` rebuilds the curved nodes FLAT and collapses Γ to the
inscribed polytope. Verified (NGSolve 6.2.2604, the order-2 hex ball):
`Curve(1)==Curve(2)==Curve(3)` all degrade the DtN floor 1.6e-5 -> 5.3e-3 (~320x)
and the Γ deviation 5.7e-5 -> 2.3e-2 (~400x), IDENTICAL for every k (the requested
order is irrelevant once the geometry is gone). `GetCurveOrder()` then reports the
new order while the geometry is flat -- a SILENT trap. (Contrast: on an OCCGeometry
mesh, which has faces to project onto, Curve(k) correctly IMPROVES Γ.) So: set the
FE order with `H1(mesh, order=k)`; change the GEOMETRY order only by re-exporting
the .vol. Guarded by the static lint rule `ngsolve-curve-after-vol-import`.
"""


def get_dtn_coarse_mesh_documentation(topic: str = "all") -> str:
    """Return DtN-spectrum / coarse-mesh-accuracy documentation."""
    topics = {
        "overview": DTN_COARSE_MESH_OVERVIEW,
        "numerics": DTN_COARSE_MESH_NUMERICS,
        "api": DTN_COARSE_MESH_API,
        "applications": DTN_COARSE_MESH_APPLICATIONS,
        "p_method": DTN_COARSE_MESH_P_METHOD,
        "formulation": DTN_COARSE_MESH_FORMULATION,
        "datasheet": DTN_COARSE_MESH_DATASHEET,
        "symmetry_hex": DTN_COARSE_MESH_SYMMETRY_HEX,
    }
    if topic == "all":
        return "\n\n".join(topics.values())
    return topics.get(topic, f"Unknown topic '{topic}'. Options: {list(topics)}")
