# TEAM 28 levitation force via 3D (axisymmetric) Cauer Ladder Network

**Result (verified 2026-06-04; external-reference validation + force-convention
fix 2026-06-20):** a 6-stage Cauer Ladder Network (CLN) reduced model
reproduces the **TEAM Problem 28 electrodynamic-levitation Lorentz force vs
height** with max |CLN - repo full-FEM| = 1.2e-6 N over the sweep.  The max
|CLN - independent reference curve| is 4.7e-4 N,
and the physically-correct levitation equilibrium (where the time-averaged
lift equals the disk weight 1.055 N) lands at **absolute disk-bottom height
z = 11.0 mm** -- matching the **published measured steady-state levitation
height z = 11.5 mm** (Karl, Fetzer, Kurz, Lehner, Rucker -- the official
TEAM 28 definition; laser triangulation, 4-measurement average) to **4%**.
TEAM 28 is a genuinely **high-Rm** levitation problem (Rm ~ 57 at the in-plane
current-loop scale R), so the lift IS the eddy reaction -- exactly the regime
where the CLN earns its keep.

**Prior art (NOT a first -- this is an open reproduction).**  The lab already
published CLN-on-TEAM-28 levitation: K. Sugahara, N. Tanimoto, Y. Takahashi,
T. Matsuo, "Cauer Ladder Network Representation with Constant Basis Functions
for Eddy Current Problems Involving Conductor Movement", COMPUMAG 2023
(Paper ID 324).  That work did MORE than this example: the full **motion-coupled
transient** levitation height `z(t)` (Matlab/Simulink, 20000 steps) compared
against the **measurement** and the conventional method, with a 4-stage CLN, in
**~7 s vs ~8 h** conventional.  It also introduced the constant-basis
`As(zgap)` expansion (`As = sum_n a_2n i_2n`) that the moving-source CLN here
and in `radia_iem_fem` builds on.  THIS example is the **open, pip-installable,
NGSolve + golden-tested reproduction** of (a frequency-domain slice of) that
result -- valuable as a reproducible open artifact, not as a first.

> **Force-convention note (the bug the published 11.5 mm caught).**  The
> reported `F_z` (and the lab `.mat` `Fz1`) is the verbatim TEAM 28 surface
> integral `Re[B_r J_t]`, which is EXACTLY **2x** the physical time-averaged
> Lorentz force `<f_z> = -(1/2) Re[J_t conj(B_r)]` (verified ratio 1.9998).
> The disk floats where the PHYSICAL lift == weight, i.e. `F_z/2 == 1.055 N`.
> An earlier version balanced the 2x integral against the 1x weight and
> reported a spurious equilibrium at dZ=+4.1 mm (absolute 14.9 mm) -- ABOVE the
> measured 11.5 mm, which is unphysical.  Comparing to the **published** height
> surfaced it; the equilibrium now uses `F_z/2` and lands at 11.0 mm.  The
> force-CONVERGENCE story (CLN vs full-FEM at a fixed height) is
> convention-independent and unaffected -- the golden still locks
> `F_z(dZ=0) = -2.1928 N`.

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

## HCurl Eddy Bubble + CLN production route

HCurl Eddy Bubble and CLN are consecutive reductions, not competing solvers.
The intended 3-D route is

    HCurl(p=6) parent
      -> face-adjacency and conductor-cycle protection
      -> EVRS spatial basis (HCurl Eddy Bubble)
      -> passive R, L, P descriptor
      -> CLN / constant-basis position interpolation
      -> force and motion coupling

The topology-aware parent reduction now covers every NGSolve HCurl volume-cell
family used by this route.  All six p=6 single-cell parents are exercised by
the implementation tests; mixed-cell meshes are inventoried per family before
the reduction.

| dimension | NGSolve cell | protected interface | reduced current | VIM interaction |
|---:|---|---|---|---|
| 3-D | TET | face | three-component `curl(T)` | affine analytic moments / exact P2 curved Duffy |
| 3-D | HEX | face | three-component `curl(T)` | six affine sub-tetrahedra |
| 3-D | PRISM / WEDGE | face | three-component `curl(T)` | three affine sub-tetrahedra |
| 3-D | PYRAMID | face | three-component `curl(T)` | two affine sub-tetrahedra |
| 2-D | TRIG | edge | out-of-plane `Jz` | planar `-log(r)/(2 pi)` interaction |
| 2-D | QUAD | edge | out-of-plane `Jz` | planar `-log(r)/(2 pi)` interaction |

`NgsolveHCurlCellFamilies` records this contract, while
`NgsolveHCurlCellVolumeInteraction` and
`NgsolveHCurlPlanarVolumeInteraction` select the dimensionally correct
epsilon-free kernel.  The 2-D path never falls back to the sampled 3-D
Laplace kernel.  Its diagnostics also report net current because an absolute
2-D inductance requires a consistent return-current/gauge convention.

The non-tetrahedral p=6 analytic path uses a tetrahedral Bernstein fit and an
exact conversion to reference monomials.  Radia's analytic Newton-potential
moments now extend through total degree 18.  The automatic family degrees are
`p-1` for TET, `2p` for WEDGE, and `3p` for HEX/PYRAMID.  The p=6 acceptance
results are:

| cell | analytic degree | projection residual | canonical sub-tets |
|---|---:|---:|---:|
| HEX | 18 | `1.90e-11` | 6 |
| WEDGE | 12 | `1.77e-13` | 3 |
| PYRAMID | 18 | `4.72e-5` | 2 |

PYRAMID modes are rational at the apex rather than finite polynomials.  The
default `1e-4` projection gate accepts the two-sub-tet degree-18 model.  A
strict `1e-8` projection gate is also available: apex-only midpoint refinement
reaches `9.34e-9` at level 8 with 114 leaf tetrahedra.  This strict setting is
an accuracy option, not the dense-Gram default.  The outer 15- and 125-point
rules on the p=6 HEX cell differ by `1.09e-4`; `outer_quad=5` is therefore the
explicit convergence check above the default 15-point rule.

Curved P2 tetrahedra no longer pass through a planarized sub-tet model.  The
curl-Piola identity is used in the form
`K(xi) = curl(T)(X(xi))*abs(det(dX/dxi))`; both physical measures are then
already contained in the two reference densities.  Radia projects `K` in the
reference tetrahedron and evaluates the Laplace distance on the exact P2 map
with the existing curved Duffy/H-matrix kernel.  On the 102-tet curved-sphere
regression at HCurl p=2, the 1020 scalar reference charges have current
projection residual `4.55e-16` and NGSolve-to-Radia geometry residual
`3.55e-16`; the reduced Gram is positive and has no diagonal epsilon.  The
production curved default uses the 125-point outer rule and eight-point
one-dimensional Duffy rule.

Warped HEX/WEDGE/PYRAMID cells retain the common analytic sub-tet kernel, but
now use uniform h refinement until both the current and geometry residuals
pass.  In the warped-HEX regression, one refinement changes 6 to 48 leaves,
reduces the current residual from `3.27e-4` to `8.63e-6`, and reduces the
geometry residual from `1.56e-2` to `6.84e-3`; the resulting Gram is positive.
The work guard is based on leaf count and the QR-compressed local polynomial
rank, so a high-degree request cannot silently expand into an impractical
H-matrix build.  The old `leaf_count^2 * monomial_count` gate is retained only
when `matrix_free=False` explicitly requests the dense verification route.

Parent orders above six are residual-controlled rather than rejected.  For a
p=7 HEX parent the formal degree is 21; the production path caps analytic
moments at degree 18 and accepts the unrefined six-sub-tet representation at a
measured residual of `3.54e-5`.  A tighter gate selects uniform h refinement.
Thus p=6 remains a studied parent order, not a hard implementation ceiling.

SIBC selection has two gates.  A conductor-air/exterior face is only a
topological SIBC *candidate*.  The half-space model is enabled only when the
local/body thickness is sufficiently larger than the skin depth.  For the
TEAM 28 aluminium disk at 50 Hz, the skin depth is 12.21 mm and the disk is
only 3 mm thick (`t/delta = 0.246`).  Therefore all 196 air-facing candidate
faces remain in the **volumetric HCurl-VIM** path; no SIBC mode is selected.

The committed structural acceptance builds a real 3-D disk mesh and obtains:

| quantity | value |
|---|---:|
| HCurl parent | p=6 |
| parent DoF | 22,814 |
| EVRS modes | 6 |
| conductor-graph cycle modes | 130 |
| selected SIBC modes | 0 |
| estimated retained modes | 136 (0.596%) |

`radia.vim.HCurlEddyCLNFromVIM` is the handoff from the spatially reduced VIM
to the passive CLN descriptor.  `radia.maglev.MovingHCurlCLNFamily` enforces
the constant-basis condition before interpolating `R(z)`, `L(z)`, and `P(z)`;
the convex interpolation preserves passivity.

The fixed-position 3-D HCurl-VIM force gate now passes without a kernel
epsilon.  Each reduced `curl(T)` mode is projected exactly to degree-5
reference polynomials, the tetrahedron self interaction uses analytic moments
from the shared degree-18 engine, and only the smooth outer integral is quadrature.  The mdx
validation produced:

| maximum tetrahedron size | p=6 parent DoF | response modes | physical `F_z` | relative error |
|---:|---:|---:|---:|---:|
| 25 mm | 22,814 | 3 | 1.101889 N | 0.513% |
| 20 mm | 28,394 | 3 | 1.098167 N | 0.173% |
| 15 mm | 40,712 | 3 | 1.092733 N | 0.322% |

The target physical force is 1.096266 N.  The largest transverse-force ratio
is 0.165%, the polynomial projection residual is below `3.4e-15`, and changing
the outer tetrahedron rule from 15 to 125 points changes the force by 0.0197%.
The structural 136-mode count is the topology-preserving a priori plan; the
3-mode count above is the final response rank after mixed-Galerkin/EVRS
compression for the three TEAM excitation ports.

The two ROM levels therefore have independent acceptance gates: the 3-D
HCurl-VIM fixed-position force passes, and the existing 25-position CLN curve
passes.  The next end-to-end gate is to drive that position sweep from the
same 3-D HCurl basis.  The generic 3-D interaction accepts affine TET, HEX,
WEDGE, and PYRAMID geometry, exact P2 curved tetrahedra, and residual-controlled
warped/curved non-tet cells.  Parent requirements above total degree 18 use an
hp degree-cap path.  Current projection, geometry, compressed scalar-charge
count, and subdivision work remain hard gates rather than silent
approximations.  See
`validation_test/maglev/team28_hcurl_vim_force_summary.json`.

HACApK is now the default end-to-end HCurl interaction operator.  On every
affine or residual-controlled TET/HEX/WEDGE/PYRAMID leaf, rank-revealing QR
selects a stable subset of the original three-component current polynomials;
it never rotates the ill-conditioned degree-18 monomial coefficients into an
artificial coefficient-space basis.  The scalar charge count is bounded by
`3 * reduced_modes * leaf_tets`, instead of `monomials * leaf_tets`.  For one
p=6 HEX this changes 7980 degree-18 monomial charges to 18 charges for one
response mode or 144 for eight modes.  The three CSR maps are registered once
in C++, and each Krylov action evaluates
`sum_c B_c^T G_HACApK B_c x` without forming the reduced inductance matrix.
Symmetrized high-order host-pair blocks cache their reverse transpose at build
time, and a self-host block reuses its one directed integration.

The exact P2 curved-tet and planar-log paths keep their existing scalar
`_ChargeGramHMatrix` and use the same composed HCurl operator, so their final
mode matrices are no longer materialized either.  `AssembleHybridVIM` preserves
these diagonal H-matrix blocks and `HybridVIMSystem.solve` uses GMRES on
`R + sL + Zs M_surface`.  Explicit `matrix_free=False`, `to_dense()`, Schur,
and mixed-Galerkin calls remain available as deliberate small-ROM verification
or condensation paths.

`HACApKSampledLaplaceInteraction` closes the former dense cross-block gap.  It
builds one stable sampled Laplace H-matrix over the volume, conductor-cycle
bridge, and surface-Omega/SIBC quadrature points, registers the three reduced
current-component maps once, and applies every reciprocal cross block through
one `sum_c B_c^T G_HACApK B_c` action.  In the default `cross_only=True`
composition, the small sampled reduced diagonal is subtracted and replaced by
the selected high-order diagonal operator.  `diagonal_bases=(volume,)` keeps an
analytic HCurl volume diagonal while the bridge and SIBC diagonals stay on the
sampled/BEM route.  This full-Gram-plus-reduced-correction construction avoids
an ACA-unstable partition-zero H-matrix and leaves no dense cross matrix in the
production Krylov path.  A projected `ngsolve.bem` backend can use the same
`build_operator(bases)` / `operator_scope` assembler contract.

This unification stops at the HCurl current space.  HDiv-MMM retains its own
BDM magnetic-charge Gram and is not inserted into the HCurl current H-matrix.
HDiv-HCurl coupling remains a separate rectangular field-coupling operator;
the two Piola maps, physical unknowns, and de Rham roles are not treated as one
isomorphic H-matrix.

## Files

| File | What |
|---|---|
| `team28_axisym_fem.py` | Repo-clean port of the lab full-FEM axisymmetric TEAM 28 solve (mixed phi-B + anisotropic-nu infinite shell). Reproduces the lab `.mat` force to **0.01%** at dZ=0. The ground-truth baseline. |
| `team28_cln_force.py`  | CLN/Cauer reduction at one height: builds K, N, F, shows the N-stage CLN force converging to full-FEM (golden). |
| `team28_cln_sweep.py`  | CLN force **vs height**, compared to the lab full-FEM `Fz1(dZ)`; recovers the physical levitation equilibrium (`F_z/2 == weight`) at absolute z ~ 11.0 mm (published 11.5 mm). |
| `validation_test/maglev/team28_hcurl_eddy_bubble.py` | Recomputes the p=6 face/cycle/SIBC policy and locks the existing 25-position full-FEM/CLN force curve as the acceptance target for the 3-D HCurl-VIM route. |
| `validation_test/maglev/team28_hcurl_vim_force.py` | Builds the p=6 3-D HCurl parent, applies topology-aware Eddy Bubble reduction, assembles the epsilon-free analytic tetrahedron VIM interaction, and verifies Lorentz force on three meshes plus an outer-quadrature check. |
| `cln_sibc_cuboid_3d.py` | Python port of the lab CLN-SIBC (Mixed Galerkin rank-(1,1) specialization) 3D cuboid core: Foster admittance + CLN reduction + Schur SIBC termination + polarizability `alpha(s)=V-Y/sigma`. The non-axisym building block. |
| `maglev_sphere_force.py` | **Isotropic induced-dipole AC levitation force** on a conducting sphere, coefficient pinned by the analytic perfect-conductor limit, frequency response reduced by CLN/Cauer. See below. |
| `ellipsoid_alpha_tensor.py` | **Shape-anisotropic polarizability tensor** of a conducting ellipsoid (analytic demag tensor + high-freq perfect-conductor `kappa_i = -V/(1-N_i)` + orientation-dependent lift). The analytic non-axisym anchor. See below. |
| `ellipsoid_alpha_omega_axisym.py` | **Full-frequency** axial `alpha_c(omega)` of a spheroid by axisymmetric FEM (uniform-field eddy solve), validated on the sphere vs `4 pi a^3 G(x)` and anchored at HF by `-V/(1-N_c)`. The eddy leg between the DC and HF analytic limits. See below. |
| `ellipsoid_alpha_tensor_3d.py` | **Transverse `m=1` tensor** via a 3D HCurl + CompactAMS solve on a graded fine-air-shell mesh. Matches the analytic sphere to ~2-3%, isotropic; the triaxial splits with the analytic ordering. Completes the tensor. See below. |
| `coil_maglev_equilibrium.py` | **Real-coil levitation equilibrium**: Radia open-boundary coil field x the verified sphere force -> stable equilibrium height + vertical stability. The Radia+NGSolve maglev workflow in action. See below. |
| `coil_sphere_eddy_force.py` | **Dipole-approximation error**: the FULL axisymmetric eddy-current Lorentz force on a sphere-in-coil vs the dipole force at `a/L ~ 0.5` (resolves the `coil_levitation` caveat). See below. |

## Isotropic levitation force (sphere) -- coefficient pinned, CLN-reduced

`maglev_sphere_force.py` builds the AC levitation force from the
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

## Shape-anisotropic ellipsoid tensor -- analytic anchor

`ellipsoid_alpha_tensor.py` makes the shape anisotropy quantitative and
verifies it analytically.  Two limits are closed-form:

- **DC**: `alpha_i(0) = 0` (no eddy response).
- **High frequency** (perfect-conductor flux exclusion): the induced
  moment along principal axis `i` is `m_i = kappa_i B_i / mu0` with
  `kappa_i = -V/(1 - N_i)`, where `N_i` is the demagnetizing factor
  (exact Osborn integral, `sum_i N_i = 1`).

| check | result |
|---|---|
| demag tensor | `sum N_i = 1`; Osborn integral == spheroid closed form (< 1e-10); sphere `1/3`; 2:1 prolate `N_c = 0.1736` |
| triaxial 5x3x1.5 mm | three DISTINCT `kappa = (-109, -130, -228) mm^3`, anisotropy `|kappa_z|/|kappa_x| = 2.09` |
| sphere limit | reduces to the verified `-2 pi a^3` exactly |
| orientation lift | `<F_z> = -(V/4mu0)/(1-N_i) grad(B0^2)` -- a body lifts **2.08x** more with `B` along its SHORT axis (most flux excluded); reduces to the verified sphere coefficient exactly |

The full-frequency eddy `alpha_i(omega)` tensor *between* `alpha=0` (DC)
and `-V/(1-N_i)` (HF) has no simple closed form for a triaxial body (the
sphere is special) -- it is computed numerically; see the axial component
below.

### Full-frequency axial `alpha_c(omega)` -- axisymmetric FEM

`ellipsoid_alpha_omega_axisym.py` computes the whole `alpha_c(omega)` curve
between those analytic anchors with an axisymmetric FEM eddy-current solve
(reusing the lab-verified TEAM 28 mixed phi-B machinery): a conducting body
of revolution in a uniform axial AC field (Dirichlet `phi = B0 r^2/2` on the
far boundary), induced moment `m_z = -pi s sigma int r phi dr dz`,
`alpha_c = m_z mu0/B0`.

| check | result |
|---|---|
| sphere validation | FEM `alpha_c(omega)` reproduces the analytic `4 pi a^3 G(x)` to **1.8%** over 2-200 kHz (BC + moment extractor correct) |
| HF anchor | reaches `-V/(1-N_c)`: sphere 1.3%, prolate 1:2 0.7%, oblate 2:1 3.0% |
| shape split | per-volume `1/(1-N_c)`: prolate 1.20 < sphere 1.48 < oblate 2.05 |
| full curve | `alpha_c` rises from DC toward the HF plateau, shape-dependent |

This is the AXIAL (`m=0`) component; the transverse (`m=1`) component
completes the full tensor (see the 3D probe below).  Convention note: the
FEM uses `s = +j omega` (engineering); `G_exact` uses `e^{-j omega t}`
(physics), so the FEM result is conjugated to match the folder's convention
(`Re[alpha]` -- the part the levitation force uses -- is unaffected).

### Transverse `m=1` component -- 3D HCurl + CompactAMS (verified)

`ellipsoid_alpha_tensor_3d.py` computes the transverse component a uniform
field along the symmetry axis (`m=0`) cannot reach: a transverse field is
`m=1`, so only a full 3D solve gives `alpha_x, alpha_y, alpha_z` at once.

| check | result |
|---|---|
| sphere magnitude | 3D `alpha_xx`, `alpha_zz` match the analytic `4 pi a^3 G(x)` to **~2-3%** (`a/delta = 1.3, 2.0`) |
| isotropy | `alpha_xx == alpha_zz` to **<0.1%** (the `m=1` transverse machinery is correct) |
| triaxial split | `5x3x1.5 mm` gives three distinct components with the analytic per-volume ordering (short axis strongest) |

**The debugging story (two BOTH-needed pieces, recorded so it is not
re-discovered):**

1. **A graded mesh with a fine AIR shell around the body.** `Re[alpha]` (the
   field-exclusion / lift part) is set by the reaction **dipole field in the
   air just outside** the conductor; a coarse air mesh under-resolves it and
   `Re` comes out ~3-4x too small (23% error), while `Im` (loss, set by the
   conductor interior) stays fine.  Refining **only the conductor does
   nothing**; refining the air shell fixes it (23% -> 2%).  This was the real
   cause -- *not* the order, *not* the formulation (an A-phi variant gave the
   identical wrong answer, because a sphere needs no scalar potential).
2. **CompactAMS + COCR** (`radia.sparsesolv_ngsolve`) to afford the resulting
   ~100k-element graded mesh -- umfpack 3D HCurl fill-in OOMs there.  order-1
   already gives 2-3% (order-2 the same), confirming air resolution, not
   order, was the lever.

Solver gotchas (recorded): build the CompactAMS preconditioner **outside**
`with TaskManager()` (nesting segfaults); CompactAMS requires `nograds=True`.

With the analytic anchors for ALL directions (`DC=0`, `HF=-V/(1-N_i)`) and
the skin-robust axial full-frequency curve above, the shape-anisotropic
polarizability **tensor is now complete and verified in every direction**.

## Real-coil levitation equilibrium -- the Radia + NGSolve workflow

`coil_maglev_equilibrium.py` composes Radia's open-boundary coil field
with the verified sphere force to find an actual levitation equilibrium --
the "Maglev = Radia + NGSolve" policy reduced to the part that needs Radia
(the sphere reaction is the analytic `G`, so no air mesh is required).  A
30 mm circular coil (`ObjFlmCur`, 10000 At) at 50 kHz lifts a 5 mm Cu
sphere; the equilibrium is the height where the lift `F_z(z) = (pi a^3/
mu0) Re[G] d|B0|^2/dz` balances gravity on the stable (descending) branch.

| check | result |
|---|---|
| Radia coil field vs analytic loop `Bz = mu0 NI R^2/(2(R^2+z^2)^1.5)` | max error 0.005% |
| sphere response | `a/delta = 16.9`, `Re[G] = -0.456` (diamagnetic lift) |
| equilibrium | `z* = 35.2 mm`, lift `46.02 mN` vs weight `46.02 mN` (residual -0.004%) |
| vertical stability | `dF/dz = -4.75 N/m < 0` (restoring; small-oscillation `f ~ 5 Hz`) |
| dipole-approximation control | `a^3` cancels in `F=weight` so `z*` depends on size only through `Re[G(a/delta)]`; shrinking the sphere drives `a/L` 0.49 -> 0.10 (clean point dipole), `z*` shifting only 15% (`Re[G]` desaturation) |

**Honest caveat** (reported in the JSON, not asserted): at `a = 5 mm` the
sphere sits where the field-gradient scale `L ~ 10 mm`, so `a/L ~ 0.5` --
the point-dipole force is approximate that close to the coil.  CHECK 5
shows the approximation is freely controllable (a smaller sphere gives a
clean `a/L << 1` at essentially the same height); a full Radia+NGSolve
eddy-current reaction-field solve would refine the force in the `a/L ~ 0.5`
regime.  That refinement is `coil_sphere_eddy_force.py` (below).

## Dipole-approximation error -- full eddy force vs the dipole force

`coil_sphere_eddy_force.py` quantifies exactly how good the point-dipole
levitation force is, by computing the FULL axisymmetric eddy-current
Lorentz force on the conducting sphere (the lab-verified TEAM 28 mixed
phi-B `(K+sN)` solver) and comparing it to the induced-dipole force at the
SAME height, using the SAME FEM coil field -- so the only difference is the
approximation.

| check | result |
|---|---|
| coil-only on-axis field vs analytic loop | 0.03% |
| `a/L = 0.50` (5 mm sphere) | `F_full/F_dipole = 1.019` (**+1.9%**: the dipole slightly under-predicts the true eddy force this close to the coil) |
| `a/L = 0.10` (1 mm sphere) | ratio `= 0.998` (**-0.2%**): the gap collapses -- confirming it is genuinely the finite-size point-dipole approximation |

**Force convention** (the load-bearing point): the full force uses the
standard time average `<f_z> = (1/2) Re[Jt conj(B_r)]`, the same
`1/2`-peak-amplitude convention as the Landau dipole coefficient, so the
two are directly comparable; the `ratio -> 1` limit as `a/L -> 0` validates
it.  The verbatim TEAM 28 integral `Re[B_r*Jt]` (a lab-`.mat` normalization)
is exactly `2x` this for the sphere (measured 2.0000; the `Im*Im` cross term
is `~6e-6`, negligible) -- a pure convention factor, not the dipole error.

## Source / provenance

- **Published benchmark (external reference)**: H. Karl, J. Fetzer, S. Kurz,
  G. Lehner, W. M. Rucker, "Description of TEAM Workshop Problem 28: An
  Electrodynamic Levitation Device", Inst. f. Theorie der Elektrotechnik,
  Univ. Stuttgart.  Official spec: Al disk R=65mm, t=3mm, m=0.107 kg; inner
  coil 960 t, outer 576 t, counter-wound; `i_hat = 20 A` peak, `f = 50 Hz`;
  rest height z=3.8mm, **measured stationary levitation height z=11.5 mm**
  (laser triangulation, 4-measurement average, Table I).  PDF in the lab
  corpus `05_TEAM_benchmark/23_problem28/`.
- Geometry + full-FEM ground truth: lab learning material
  `W:\00_CAE\NGSolve\01_菅原\2024_08_TEAM28` (axisymmetric NGSolve TEAM 28,
  DC / 50Hz / 50Hz_可動 / Transient + field-validation figure).
  Disk: Al, R=65mm, t=3mm, sigma=3.4e7; coils: 960t/+20A (r=41mm) and
  576t/-20A (r=87.5mm) counter-wound, 50 Hz.
- CLN theory: `radia_mcp.mor` (mor_cln); Kameari-Ebrahimi-Sugahara-
  Shindo-Matsuo 2018, IEEE TMag 54(3):7201804.
- **Prior CLN-on-TEAM-28 (the result this reproduces)**: K. Sugahara,
  N. Tanimoto, Y. Takahashi, T. Matsuo, "Cauer Ladder Network Representation
  with Constant Basis Functions for Eddy Current Problems Involving Conductor
  Movement", COMPUMAG 2023 (Paper ID 324) -- motion-coupled transient
  levitation height vs measurement + conventional method, 4-stage CLN, ~7 s
  vs ~8 h.
- Method context: `radia_mcp.maglev` topics `cln_mor_control` /
  `radia_iem_fem`; the CLAUDE.md policy "Maglev Analysis: Radia + NGSolve".

## Run

```bash
python team28_axisym_fem.py      # full-FEM baseline  -> -2.1925 N @ dZ=0
python team28_cln_force.py       # CLN convergence    -> 5-stage golden
python team28_cln_sweep.py       # CLN force vs height -> physical equilib z~11.0mm (pub 11.5mm)
python validation_test/maglev/team28_hcurl_eddy_bubble.py  # run from repo root
python cln_sibc_cuboid_3d.py     # CLN-SIBC 3D cuboid core (alpha, Schur-F)
python maglev_sphere_force.py  # isotropic levitation force, coeff pinned
```
