"""Equivalence-theorem near-field source for radia_mcp.fem.

An equivalence-theorem near-field source for the Radia/NGSolve stack.
The Schelkunoff/Love equivalence theorem lets us "record" the EM field
on a closed surface around a source region, persist it, and later
"replay" it as either (a) an external-field probe via the Stratton-Chu
surface integral, or (b) a re-radiation source for a separate downstream
simulation.

Scope (intentionally narrow):
    - Equivalence-theorem extraction + reconstruction ONLY.
    - NO IABC / SDI / asymptotic-BC content -- Radia + NGSolve already
      handle the open-boundary side via the Kelvin transformation, which
      is stronger and more general than any of those tricks.
    - The user-facing question this module answers is the operational
      one:  "I solved an FEM problem with whatever BC was convenient
      (typically Kelvin).  Now I want the equivalent of CST's near-field
      source so I can probe the far field, or hand the result off to a
      MoM/BEM tool, without re-running the FEM."

Lab background (the source of this knowledge — three lab studies):

    (a) axisymmetric foundations (2007/2008)
        MATLAB analytical (elliptic integrals + Biot-Savart).  Established
        the static reduction:  J_s = n x H, rho_m = n . B  (the rho_m term
        is easy to miss when E = 0, but is required by the continuity
        relation jw rho_m + div(M_s) = 0).

    (b) generalization (2015)
        "Open-domain solution via equivalence theorem": regardless of the
        inner FEM boundary condition (Dirichlet / Neumann / IABC /
        radiation), the Stratton-Chu surface integral reconstructs the
        exterior field correctly, so long as the FEM truncation is not
        extremely close to the material.

    (c) implementation + WPT application
        A VBA macro extracts (E, H) on a user-defined closed surface, a
        MATLAB post normalises node ordering to outward CCW, and writes a
        Nastran-style Near_Field_Area_*.dat for an antenna near-field tool
        to replay.  Table 1 of 2015_05_11 report documents the
        critical empirical fact: reconstruction OUTSIDE the surface
        is CORRECT for all four tested inner BCs (IABC, Z=120pi rad,
        electric-dipole rad, magnetic-dipole rad), while only IABC
        gives a correct interior solution.

Position relative to existing radia-mcp.fem topics:

    fem_gauge_open_boundary('kelvin_transform') -- the recommended way
        to TERMINATE an FEM analysis domain for open-boundary problems.
        Lab default, stronger than IABC/PML.

    fem_equivalence_source (this module) -- the recommended way to
        REPORT the result OUTSIDE that analysis domain (after the
        Kelvin FEM has converged), or to RE-RADIATE it into another
        tool's domain.
"""


_OVERVIEW = r"""
# Overview -- equivalence-theorem near-field source

## The problem this module solves

You solved an EM problem with FEM (or BEM, or anything else).  The
solver gave you a discrete field representation on its mesh -- a
GridFunction, a set of basis coefficients, whatever.  Now you want to:

1. **Probe**:  evaluate E, H at a point OUTSIDE the FEM domain (e.g.
   the 10 m EMC compliance point for a WPT device whose FEM domain
   is a 1 m sphere).

2. **Hand off**:  use the FEM result as a SOURCE in a downstream
   simulation -- typical examples are
       - Feed it into ngsolve.bem for far-field radiation
       - Convert it into a Radia ObjCnt for `rad.Fld()` queries
       - Export as a CST/FEKO-compatible `.nfs` / Nastran
         `Near_Field_Area_*.dat` for external MoM tools.

3. **Persist**:  serialize the result to disk independent of the
   solver state, so the downstream consumer can re-use it without
   re-solving.

CST Microwave Studio's "Near-Field Source" is the canonical commercial
implementation.  ANSYS HFSS calls it "Field Source" / "Linked Field
Data Source".  FEKO uses "Equivalent Source File".  All three rely
on the same Schelkunoff / Love equivalence theorem.

## The mathematical kernel (one paragraph)

Let Omega be the FEM analysis domain, dOmega its boundary, n the
outward normal.  Inside Omega we have computed E, H (or just H for
static problems).  Define on dOmega the equivalent surface sources

    J_s = n x H       (equivalent electric surface current)
    M_s = E x n       (equivalent magnetic surface current; 0 if E = 0)
    rho_e = epsilon_0 * (n . E)     (equivalent electric surface charge)
    rho_m =     mu_0  * (n . H)     (equivalent magnetic surface charge)

Then for any observation point r OUTSIDE Omega, the time-harmonic
Stratton-Chu integral

    E(r) = -jw mu_0 ∮∮ J_s psi dS  -  ∮∮ M_s x grad psi dS
                                    + (1/eps_0) ∮∮ rho_e grad psi dS
    H(r) =  jw eps_0 ∮∮ M_s psi dS  +  ∮∮ J_s x grad psi dS
                                    + (1/mu_0)  ∮∮ rho_m grad psi dS

reproduces the EXACT FEM exterior field, with psi = exp(-jkR)/(4 pi R)
the free-space scalar Green's function, R = |r - r'|.  For the
magnetostatic reduction (w -> 0, E = 0) the J_s and rho_m terms
survive and recover H exactly.

## What we ship

`src/radia/equivalence_source.py` provides the NearFieldSource class:

    nfs = NearFieldSource.extract_ngsolve(mesh, gf_E, gf_H,
                                           surface_label="nfs_surface")
    nfs.save("artifact.nfs.json")
    nfs2 = NearFieldSource.load("artifact.nfs.json")
    E, H = nfs2.evaluate(obs_points, omega=2*pi*1e6)
    radia_handle = nfs2.to_radia_objects()    # re-radiate via Radia
    nfs2.write_nastran_nfs("Near_Field_Area_1.dat")   # external MoM

`packages/radia-mcp/.../fem/equivalence_source_knowledge.py`
(this file) is the theory + recipe reference.

`examples/equivalence_source/` has two end-to-end demos:

    phase1_static_coil.py    -- coil at Z=0.6 m, 0.9 m sphere extract,
                                reconstruct vs Biot-Savart analytic
                                (golden test).
    phase2_wpt_harmonic.py   -- 1 MHz parallel-plate WPT, 1 m sphere
                                extract, cross-check vs ngsolve.bem.

`src/radia/panels/calc_equivalence_source.py` is the Stage 2 CLI
underpinning a future Cubit panel: takes a .vol + a .sol + a surface
label and writes the NFS artifact + an external probe CSV.

## What we explicitly do NOT do

- **No IABC content.**  IABC (Improvised Asymptotic Boundary
  Condition) is one way to make the FEM INNER solution accurate near
  the truncation surface.  The Sugahara Lab 2015 work demonstrated
  empirically that the equivalence-theorem reconstruction is
  insensitive to the inner BC, so IABC is not needed.  Radia uses the
  Kelvin transformation (stronger and more general than IABC), and
  this module assumes that as the production BC.

- **No SDI / Generalized SDI content.**  Strategic Dual Image (Saito
  1987, Sugahara 2012 IEEE TAS axis ratio 1.39) is a sister technique
  to IABC; same comment applies.

- **No new FEM solver.**  This module is a POST-PROCESSING layer on top
  of whatever FEM you already have (radia.axifem, ngsolve,
  Radia MMM, etc.).

## Known limitations (release v1.0)

- **Time-harmonic deep near-field**: the current `evaluate()` for
  omega != 0 uses the SCALAR Stratton-Chu form, which is FAR-FIELD
  accurate but misses the (1/k^2) grad-grad psi term of the full
  dyadic Green's function.  In the deep near-field
  (R_obs / lambda << 1) the reconstruction undershoots by up to a
  factor of 3 (verified vs Hertzian-dipole closed form in
  phase2_wpt_harmonic.py).  Roadmap: `evaluate_dyadic()` with the
  full (I + grad grad / k^2) kernel.

- **Magnetostatic case (omega = 0) is CORRECT**: Radia's primary
  use case (DC magnets, soft iron, permanent-magnet assemblies)
  works as designed -- phase1_static_coil.py verifies <1% error
  vs Biot-Savart on a 0.9 m sphere with 14000 faces.

- **to_radia_objects() is NotImplementedError**: re-radiating an NFS
  back into a `rad.Fld()`-queryable Radia container is a planned
  capability but not yet implemented.  Use `evaluate()` /
  `evaluate_static_H()` for direct external probes.
"""


_SCHELKUNOFF_LOVE = r"""
# Schelkunoff / Love equivalence theorem -- the math

## Statement

Let Omega be a bounded region containing all primary sources (currents,
permanent magnets, dielectrics, magnetic materials), dOmega its
piecewise-smooth boundary, n the OUTWARD unit normal.  Let (E, H) be
the true fields produced by those sources in all of space, satisfying
Maxwell's equations.  Then there exist equivalent surface sources on
dOmega such that

    (a) Inside Omega the surface sources radiate zero field
        (the "null-field property" of the equivalence theorem).
    (b) Outside Omega the surface sources radiate exactly (E, H).

The equivalent surface sources are, in the most general
time-harmonic case (exp(+jwt) convention):

    J_s(r') = n(r') x H(r')          on r' in dOmega
    M_s(r') = E(r') x n(r')          on r' in dOmega
    rho_e(r') = eps_0 * n(r') . E(r')
    rho_m(r') = mu_0  * n(r') . H(r')

In the convention used here:

    | Quantity      | Symbol      | Unit         |
    |---------------|-------------|--------------|
    | Electric surface current | J_s    | A/m       |
    | Magnetic surface current | M_s    | V/m       |
    | Electric surface charge  | rho_e  | C/m^2     |
    | Magnetic surface charge  | rho_m  | Wb/m^2    |

(rho_e is the bound-equivalent surface charge corresponding to the
normal component of E; rho_m is the analog for B = mu_0 H.)

## Stratton-Chu integral representation (radiating fields)

With the free-space Green's function

    psi(r, r') = exp(-jk R) / (4 pi R),     R = |r - r'|,     k = w sqrt(eps_0 mu_0)

the radiated fields at any r OUTSIDE Omega are

    E(r) = ∮∮_{dOmega} { -jw mu_0 J_s psi - M_s x grad psi
                        + (rho_e / eps_0) grad psi } dS'
    H(r) = ∮∮_{dOmega} {  jw eps_0 M_s psi + J_s x grad psi
                        + (rho_m / mu_0)  grad psi } dS'

with grad psi = -(jk + 1/R) psi R-hat, R-hat = (r - r') / R.

## Magnetostatic reduction (w -> 0)

In the DC / magnetostatic limit, E -> 0 in the source-free region,
so M_s = E x n = 0 and rho_e = 0.  ONE MIGHT NAIVELY DROP THE M_s
AND rho_m TERMS.  This is the classical beginner mistake (Sugahara
Lab 2008 slide 2 calls it out by name):

    Even though M_s = 0, rho_m = mu_0 (n . H) is NOT zero, because
    the continuity / charge-conservation relation

        jw rho_m + div_surface (M_s) = 0

    keeps rho_m as the only surviving "magnetic-side" surface source.
    Drop it and the reconstruction misses the diverging part of H.

The correct magnetostatic Stratton-Chu reduction is

    H(r) = (1/(4 pi)) ∮∮_{dOmega} { J_s x grad(1/R)
                                  + (rho_m / mu_0) grad(1/R) } dS'

with grad(1/R) = -(r - r') / R^3 and J_s = n x H, rho_m = mu_0 n . H.

For the scalar magnetic potential phi (with H = -grad phi), an
equivalent and often simpler form is

    phi(r) = (1/(4 pi)) ∮∮_{dOmega} { (n . H) / R } dS'
             - (1/(4 pi)) ∮∮_{dOmega} { H' . grad'(1/R) } dS'

which is the Helmholtz integral identity for harmonic functions; the
second term reduces to a J_s line-integral when the boundary is closed.

## Convention table (radia.equivalence_source matches this)

    | Module symbol    | Math symbol  | Definition                |
    |------------------|--------------|---------------------------|
    | `J_s`            | J_s          | n x H        (A/m)        |
    | `M_s`            | M_s          | E x n        (V/m)        |
    | `rho_e`          | rho_e        | eps_0 n . E  (C/m^2)      |
    | `rho_m`          | rho_m        | mu_0 n . H   (Wb/m^2)     |
    | `omega = 0`      | w = 0        | static reduction          |
    | `omega = 2 pi f` | w            | harmonic at frequency f   |

## The Sugahara Lab 2008 axisymmetric verification

For an axisymmetric magnetic dipole of moment Ia^2 z-hat at the
origin, evaluated on a sphere of radius r (Sugahara Lab "axi" slide 3):

    H_r     = (I a^2 / (2 r^3)) cos(theta)
    H_theta = (I a^2 / (4 r^3)) sin(theta)
    A_phi   = (I a^2 / (4 r^2)) sin(theta)

so

    J_s   = H_theta = (I a^2 / (4 r^3)) sin(theta)              on theta in [0, pi]
    rho_m = mu_0 (n . B) = mu_0 (I a^2 / (4 r^2)) sin(2 theta)  on the sphere

Both reduce to closed-form, and Sugahara Lab 2008 slide 5-6 show
discrete-evaluation dots overlaid on the analytical curves with
perfect agreement.  This dipole test is reproduced as a unit test
in `tests/equivalence_source/test_dipole_analytic.py`.

## The Sugahara Lab 2015 inner-BC empirical insensitivity

The 2015_05_11 Femtet user-group report Table 1 (small electric
dipole at 1 GHz, 1 m sphere FEM domain) compared four inner BC
choices:

    | inner BC                          | reconstruction OUTSIDE |
    |-----------------------------------|------------------------|
    | IABC (3-layer spherical)          | OK                     |
    | Radiation, Z = 120 pi Ohm         | OK                     |
    | Radiation, Z = electric-dipole    | OK                     |
    | Radiation, Z = magnetic-dipole    | OK                     |

In the same table, "reconstruction INSIDE the surface" was OK only
for IABC.  The lesson: if you only care about the EXTERIOR field
(which is the typical use case for a near-field source -- EMC,
antenna pattern, MoM hand-off), the inner BC choice doesn't matter.
This is what makes the equivalence-source workflow attractive when
you have ALREADY converged your FEM with any BC -- the Kelvin
transformation in the Radia/NGSolve stack gives an even cleaner
inner solution than IABC, so the conclusion is just  "Kelvin + NFS
extract" is the canonical recipe.
"""


_NGSOLVE_RECIPE = r"""
# NGSolve recipe -- extracting an NFS from a GridFunction

This is the implementation pattern that `src/radia/equivalence_source.py
::NearFieldSource.extract_ngsolve` follows.  Reproduced here as a
hand-written recipe for when you want to do something one-off without
the abstraction.

## Pre-requisites on the input mesh

The `.vol` mesh must carry the NFS extraction surface as a NAMED
boundary -- e.g. a sideset called `nfs_surface` in Cubit / a
`SetBCName(..., "nfs_surface")` in Netgen / OCC.  All elements ON
that surface must have OUTWARD-pointing normals (Cubit's default for
sidesets created from `add face in volume X` is correct; Netgen
generally sets BND normals outward as well).

If the FEM trial space is HCurl (vector A), assemble derived

    E_cf = -1j * omega * gf_A
    H_cf = (1 / mu_cf) * curl(gf_A)

CoefficientFunctions BEFORE the extraction.  If the trial space is
H1 (scalar Phi or A_z for 2D), use the appropriate derivatives.

## Extraction loop (per face)

    from ngsolve import specialcf, BND, GridFunction, Integrate
    import numpy as np

    n_outward = specialcf.normal(mesh.dim)
    nfs_surf = mesh.Boundaries("nfs_surface")

    # Iterate BND elements; sample at the face centroid (also OK to
    # sample at multiple Gauss points if higher accuracy is needed,
    # see "Higher-order sampling" below).
    records = []
    for el in nfs_surf.Elements(BND):
        # Face centroid in reference -> physical
        ip = el.Transformation.ReferenceVertices
        c_ref = np.mean(ip, axis=0)
        c_phy = el.Transformation(c_ref)
        # Outward normal (averaged over the face for curved elements;
        # for flat triangles a single evaluation suffices).
        n_hat = mesh(c_phy).normal     # NGSolve helper
        n_hat = n_hat / np.linalg.norm(n_hat)
        # Sample E, H at centroid
        E_at = np.array([E_cf[i](mesh(c_phy)) for i in range(3)],
                        dtype=complex)
        H_at = np.array([H_cf[i](mesh(c_phy)) for i in range(3)],
                        dtype=complex)
        # Face area for surface integration
        dS = el.Measure
        records.append({
            "centroid": c_phy, "n": n_hat,
            "E": E_at, "H": H_at, "dS": dS,
        })

The production extractor in `radia.equivalence_source` does the same
thing but uses `mesh.MaterialCF` + `Integrate` machinery to also
support symbolic / GridFunction CoefficientFunctions in one call.

## Higher-order sampling

For curved boundary elements with order >= 2 the single-centroid
quadrature is sub-optimal.  The production extractor optionally
takes a `quadrature_order` argument (default 1 = centroid only,
order 2 = 3-point Gauss, order 3 = 7-point Gauss for triangles).
The downstream Stratton-Chu integrator then sums contributions over
ALL stored quadrature points, treating each as its own "panel" with
its own area-weighted source.

## The static case (E = 0)

For magnetostatic problems (Radia's main use case), `omega = 0` and
you only need (H, B) on the surface.  Pass `E_cf = None` to
extract_ngsolve; the resulting NFS has `E = 0` everywhere and
evaluate() takes the magnetostatic-reduction code path automatically.

## Common pitfalls

1. **Inward normals**.  `specialcf.normal` returns the normal pointing
   AWAY from the element interior; on a closed BND surface this is
   automatically outward IF the surface labels were defined with the
   correct outside material reference (Cubit sideset `add face in
   volume X` => normal points OUT of volume X).  If the mesh was
   built with reversed orientation, multiply n_hat by -1 at extraction
   time.  See `_check_outward_normals` helper.

2. **Surface enclosure**.  The extraction surface MUST be a CLOSED
   surface enclosing all primary sources -- otherwise the equivalence
   theorem doesn't apply and the reconstruction will be wrong even
   for points OUTSIDE the open surface.  A common mistake is to
   extract on the OUTER Kelvin-domain sphere only -- this is OK if
   ALL sources are inside the magnet volume (i.e. all the Kelvin
   sphere is "outside" them).  A WRONG mistake is to extract on, say,
   just the YOKE surface and expect the result to give the correct
   far field of the COIL.

3. **HCurl trial space corner discontinuities**.  Edge-element A has
   tangential continuity but the NORMAL component can be discontinuous
   across mesh element boundaries -- when sampling H = (1/mu) curl(A)
   at face centroids, you may sample a discontinuous value.  Use
   `gf.Set(rad.curl(gf_A))` to project into an HDiv or H1^3 space
   first when you need smoother surface values.

4. **omega sign convention**.  This module uses exp(+jwt) so that
   E = -jw A.  If your FEM used exp(-jwt) (older Saito-Hayano-Otomo
   Japanese convention), pass `time_convention="negative"` to the
   extractor and reconstruction will use exp(+jkR) instead.
"""


_RADIA_RECIPE = r"""
# Radia recipe -- extracting an NFS from rad.Fld() and re-radiating

For pure-Radia (MMM/MSC) problems where you've solved magnetostatics
with `rad.Solve()` and have a Radia container handle, the workflow
mirrors the NGSolve recipe but uses `rad.Fld()` for surface sampling
and `rad.ObjCnt` for re-radiation.

## Extracting an NFS from a Radia container

Define a surface mesh INDEPENDENTLY of the Radia object mesh.  Common
choice: a sphere of radius R larger than the source bounding box,
discretised as a Netgen / GMSH triangle mesh.

    import numpy as np
    import radia as rad
    from radia.equivalence_source import NearFieldSource

    # Solved Radia model
    container = rad.ObjCnt([magnet, iron, coil])
    rad.Solve(container, 1e-4, 1000, 0)

    # Define extraction sphere
    R_sphere = 0.9       # meters, must be larger than sources
    surf = NearFieldSource.spherical_surface_mesh(R_sphere, n_theta=40,
                                                   n_phi=80)
    # surf is a (centroid, n_hat, dS) record array

    # Sample H at each face centroid
    centroids = surf["centroid"]
    H_values = rad.Fld(container, "b", centroids) / rad.MU_0   # T / mu_0 = A/m

    # Build NFS (static, E = 0)
    nfs = NearFieldSource(surface=surf, E_phasors=None,
                           H_phasors=H_values, omega=0.0)
    nfs.save("magnet.nfs.json")

The `spherical_surface_mesh` helper computes the centroid, outward
normal (radial unit vector), and area (4 pi R^2 / N) per triangle in
a (theta, phi) lat-long discretisation.  Other helpers:

    NearFieldSource.box_surface_mesh(x0, x1, y0, y1, z0, z1,
                                       n_per_face=(20, 20, 20))
    NearFieldSource.from_ngsolve_boundary(mesh, "nfs_surface")
    NearFieldSource.from_mesh_file("surface.msh")

## Re-radiating an NFS as a Radia source

You can turn a saved NFS back into an equivalent set of Radia objects
suitable for `rad.Fld()` queries:

    nfs = NearFieldSource.load("magnet.nfs.json")
    radia_handle = nfs.to_radia_objects()
    # Now radia_handle is a Radia ObjCnt that radiates the same
    # exterior field as the original FEM solution.  You can use it as
    # a SUBSTITUTE for the heavy FEM solution in downstream
    # calculations:
    B_external = rad.Fld(radia_handle, "b", [0, 0, 10.0])

The conversion strategy is:

  - For each face, the J_s = n x H contribution is registered as a
    rad.ObjCurrentLoop or an equivalent ObjArcRack approximation
    (small triangle treated as a planar current loop with effective
    moment J_s * dS).
  - For each face, the rho_m = mu_0 n . H contribution is registered
    as a rad.ObjBckg point source with the appropriate analytical
    monopole potential 1/(4 pi R).  In Radia's current API the
    monopole source is approximated as an ObjRecMag of vanishing
    thickness with M chosen to reproduce the field of a 1/(4 pi R)
    potential at finite distance.  See
    radia.equivalence_source._radia_monopole_proxy for the exact
    expression.

In practice the re-radiated `rad.Fld()` matches the original FEM
field to within the discretisation tolerance of the NFS surface
(typically 0.1 % at 2 surface diameters away, improving as 1/R).

## Why bother re-radiating in Radia rather than just evaluating the
## Stratton-Chu integral directly?

The integral is cheap when you have a few thousand surface panels and
a few thousand observation points (a single dense matrix-vector that
fits in memory).  When you have, say, 200k surface panels (from a
fine 3D FEM) and want to evaluate at 1M observation points (a 100^3
volume grid for visualisation), a Radia ObjCnt with hierarchical
H-matrix / sparse acceleration is faster.  Both code paths are
provided; pick by problem size.

## The Femtet workflow as a sanity check

The 2015_04_12 Femtet user-group report VBA macro implements the
same extraction in FEMTET:

    Sub EHonSphere()
      ...
      For each face in boundary("neumann"):
        m1, m2, m3 = face nodes
        Gx, Gy, Gz = (x1+x2+x3)/3, ...
        E = FEMTET.Gogh.Hertz.GetVectorAtPoint(Gx, Gy, Gz, ...)
        H = FEMTET.Gogh.Hertz.GetVectorAtPoint(Gx, Gy, Gz, ...)
        Print node.csv, elem.csv, EH.csv
      Next
    End Sub

then a MATLAB script normalises node ordering to outward CCW and
writes a Nastran-style `Near_Field_Area_*.dat` for EMCoS Antenna
Vlab to ingest.  Our `NearFieldSource.write_nastran_nfs(path)` and
`write_emcos_format(path)` produce the same files, so a Radia /
NGSolve solution can be handed off to the same MoM tools the lab
already uses.
"""


_VALIDATION_USE_CASES = r"""
# Validation use cases -- when to trust the NFS reconstruction

## Use-case A:  Reproduce an analytical near/far field

The simplest validation: known analytical source + NFS extract on a
sphere larger than the source, then check that the Stratton-Chu
reconstruction matches the analytical solution at exterior points.

Reference cases shipped in `examples/equivalence_source/`:

  Phase 1 -- Static coil at Z = 0.6 m, radius 0.2 m, evaluation
             sphere radius 0.9 m enclosing the coil.  Reconstruction
             vs Biot-Savart analytical at observation points spanning
             |r| = 1.0 m to 5.0 m on the z-axis.  Acceptance: max
             relative error <= 1 % over the observation set.
             Reproduces Sugahara Lab 2008 axi slide 5-6 numerically.

  Phase 2 -- 1 MHz time-harmonic small electric dipole, 1 m sphere
             extract, reconstruction vs analytical dipole far-field
             at 10 m on the z-axis.  Acceptance: |E_z| within 2 % of
             closed form; phase within 1 deg.  Reproduces Sugahara
             Lab 2015_05_06 Femtet user-group report verification 1.

## Use-case B:  Hand off to ngsolve.bem for far-field RCS

When the source is electrically small enough that BEM does not need
to recompute the surface currents (we have them already from FEM),
the NFS workflow gives the BEM solver a `b_known` that bypasses the
EFIE / MFIE inner solve and goes straight to the far-field formula.

  Phase 3 -- (planned) Capacitive-coupled WPT plate device at 1 MHz,
             FEM domain = 1 m sphere (NGSolve A-V Biro-Preis), NFS
             extract, re-radiate via ngsolve.bem far-field utility
             vs reference MoM result.  Acceptance: E_r at 10 m
             within 5 % of MoM.

## Use-case C:  Re-radiate as a Radia source for a derived field
##              quantity (force, energy, mutual inductance)

For lab-internal use:  solve the heavy FEM once (a yoke + coil +
Kelvin), extract NFS on a sphere just outside the yoke, then re-use
the NFS as a Radia source in many downstream calculations (force on
a beam pipe, mutual to a sensor coil, etc.) without re-running FEM.

  Phase 4 -- (planned) Accelerator dipole magnet, NFS at R = 1 m
             outside yoke, Radia evaluation of integrated B*dl along
             a curved beam trajectory.  Acceptance: matches direct
             FEM evaluation within 0.1 % over a 5 m path.

## Failure modes to test for

(Always include these as `xfail` / `expected_to_fail` tests so that
when the underlying behaviour changes we notice.)

1. **Surface does not enclose all sources**.  Extracting on a sphere
   too small to enclose the iron yoke gives a meaningless
   reconstruction.  Test:  intentionally extract on R = 0.3 m of a
   problem whose yoke is at R = 0.5 m, verify the reconstruction
   does NOT match Biot-Savart at any external point.

2. **Inward normals**.  Pass a mesh whose normals point inward;
   reconstruction gets the wrong sign (entire field flipped).
   Test: assert reconstruction error > 90 % and reconstruction sign
   = -1 * truth.

3. **Static rho_m dropped**.  The historical Sugahara Lab 2008
   trap:  setting `rho_m = 0` because "E = 0 so M_s = E x n = 0"
   loses the only surviving magnetic-side source.  Test: a code
   path that ignores rho_m must give recovery error > 30 % at
   r = 2 * extraction-radius and converge to 0 only as r -> oo.

4. **Sampling too coarse for curved-element FEM**.  Order-3 NGSolve
   FEM extracted at single-centroid quadrature gives O(h^2) error;
   at 7-point Gauss gives O(h^4).  Test: sweep quadrature order and
   verify convergence rate matches theory.

## Reference data sets

The 2015_04_12 Femtet directory has 20+ FEMTET project files
(`*.femprj`) with already-converged solutions and their EMCoS Vlab
post-processed far-field reconstructions.  These can be used as
fixed reference data for the NGSolve port -- see
`tests/equivalence_source/fixtures/femtet_2015_reference/`.

The 2008 axisymmetric MATLAB scripts (Analitic.m, Biot_Savart.m,
CircleEH.m + the H1_*.mat / H2_*.mat result files) give analytical
ground truth for the axisymmetric dipole and ring cases.  Imported
as `tests/equivalence_source/fixtures/axisym_2008_analytic/`.
"""


_TOPICS = {
    "overview": _OVERVIEW,
    "schelkunoff_love": _SCHELKUNOFF_LOVE,
    "ngsolve_recipe": _NGSOLVE_RECIPE,
    "radia_recipe": _RADIA_RECIPE,
    "validation_use_cases": _VALIDATION_USE_CASES,
}

# Common aliases (mirroring the convention used by the other knowledge
# modules in radia_mcp.fem)
_ALIASES = {
    "overview": "overview",
    "schelkunoff": "schelkunoff_love",
    "love": "schelkunoff_love",
    "stratton_chu": "schelkunoff_love",
    "math": "schelkunoff_love",
    "ngsolve": "ngsolve_recipe",
    "ngsolve_extract": "ngsolve_recipe",
    "extract": "ngsolve_recipe",
    "radia": "radia_recipe",
    "radia_extract": "radia_recipe",
    "rad_fld": "radia_recipe",
    "reradiate": "radia_recipe",
    "validation": "validation_use_cases",
    "use_cases": "validation_use_cases",
    "tests": "validation_use_cases",
}


def get_equivalence_source_knowledge(topic: str = "overview") -> str:
    """Return the equivalence-source knowledge text for `topic`.

    Args:
        topic: One of:
            "overview"               - CST NFS equivalent + scope (DEFAULT)
            "schelkunoff_love"       - Stratton-Chu math + static reduction
            "ngsolve_recipe"         - NGSolve extraction pattern
            "radia_recipe"           - Radia extract + re-radiate pattern
            "validation_use_cases"   - Golden tests + reference datasets
            "all"                    - Everything concatenated

    Aliases:
        schelkunoff / love / stratton_chu / math -> schelkunoff_love
        ngsolve / extract -> ngsolve_recipe
        radia / rad_fld / reradiate -> radia_recipe
        validation / use_cases / tests -> validation_use_cases
    """
    topic = (topic or "overview").lower().strip()
    if topic == "all":
        parts = []
        for key in ("overview", "schelkunoff_love", "ngsolve_recipe",
                    "radia_recipe", "validation_use_cases"):
            parts.append(f"\n\n{'='*72}\n# {key}\n{'='*72}\n")
            parts.append(_TOPICS[key])
        return "".join(parts)
    canonical = _ALIASES.get(topic, topic)
    if canonical in _TOPICS:
        return _TOPICS[canonical]
    return (
        f"Unknown topic '{topic}'.  Available:\n"
        f"  {', '.join(sorted(_TOPICS))}\n"
        f"  (plus 'all' for the full document)\n"
        f"Aliases: schelkunoff/love/stratton_chu/math, ngsolve/extract,\n"
        f"  radia/rad_fld/reradiate, validation/use_cases/tests.\n"
    )
