"""XFEM in COMSOL Multiphysics: theory, workflow, and worked examples.

Knowledge layer for the eXtended Finite Element Method (XFEM) based on
the lab reference paper:

    Ahmad Jafari, Pooyan Broumand, Mohammad Vahab, Nasser Khalili,
    "An eXtended Finite Element Method Implementation in COMSOL
    Multiphysics: Solid Mechanics," arXiv:2109.03153v1 [cs.CE], Sep 2021.
    School of Civil and Environmental Engineering, UNSW Sydney.
    Companion code: https://github.com/ahmadjafari93/xfem-comsol

This module covers:
    - XFEM theory: partition of unity, Heaviside + asymptotic
      crack-tip enrichment, level-set crack tracking
    - COMSOL-specific implementation pattern: 6 Solid Mechanics
      modules (1 standard + 1 Heaviside + 4 tip) coupled through
      modified stress / strain definitions in the Equation View
    - Mixed-mode crack growth: K_I / K_II via interaction integral,
      maximum hoop stress propagation criterion
    - Worked examples: center crack (verification), randomly
      distributed cracks, plate with hole, double-hole multi-crack,
      penny-shaped 3D crack

Lab paths to the source files are exposed via LAB_FILES so that
downstream automation (e.g. subprocess COMSOL batch) can locate the
.mph models without re-hardcoding.

All text is ASCII (cp932-safe) and in English.
"""

from __future__ import annotations


# ---------------------------------------------------------------------------
# Reference paper file paths (lab-internal LAN share)
# ---------------------------------------------------------------------------

LAB_FILES = {
    "root":             r"S:/COMSOL/2022_05_12_XFEM/",
    "paper_pdf":        r"S:/COMSOL/2022_05_12_XFEM/An eXtended Finite Element Method Implementation in COMSOL Multiphysics Solid Mechanics.pdf",
    "readme":           r"S:/COMSOL/2022_05_12_XFEM/README.md",
    "license":          r"S:/COMSOL/2022_05_12_XFEM/LICENSE",
    # The two .mph examples are large; do NOT open as text.
    "penny_3d_mph":     r"S:/COMSOL/2022_05_12_XFEM/XFEM 3D penny-shaped crack.mph",
    "mixed_mode_mph":   r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/XFEM mixed-mode crack growth in plate with hole.mph",
    "mixed_mode_dir":   r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/",
    "mixed_mode_archive_zip": r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole.zip",
    "infinite_archive_rar":   r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack infinite domain.rar",
    # MATLAB support scripts that drive level-set updates from COMSOL
    "preprocess_m":     r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/preprocess.m",
    "phi_m":            r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/phi.m",
    "interpol_m":       r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/interpol.m",
    "update_crack_m":   r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/update_crack.m",
    "read_crack_m":     r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/read_crack.m",
    "last_angle_m":     r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/last_angle.m",
    "p_poly_dist_m":    r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/p_poly_dist.m",
    "mesh_mphtxt":      r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/mesh.mphtxt",
    "read_me_txt":      r"S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/READ ME.txt",
}


# ---------------------------------------------------------------------------
# Module-level knowledge blocks
# ---------------------------------------------------------------------------

OVERVIEW = r"""
# XFEM in COMSOL Multiphysics: overview

Reference paper (lab archive):
    Jafari, Broumand, Vahab, Khalili, 2021
    "An eXtended Finite Element Method Implementation in COMSOL
    Multiphysics: Solid Mechanics", arXiv:2109.03153v1 [cs.CE], Sep 2021.
    (UNSW Sydney; companion code https://github.com/ahmadjafari93/xfem-comsol)

XFEM (eXtended FEM) enriches a standard finite element approximation
space with extra functions that capture features the mesh cannot
naturally resolve:

  - jumps across cracks / material interfaces (Heaviside enrichment),
  - 1/sqrt(r) stress singularities at crack tips (asymptotic
    enrichment using the four Westergaard near-tip functions),
  - other features built into the local solution (voids, slip
    planes, dislocations, weak discontinuities).

The mathematical scaffolding is the Partition of Unity (PU) method
(Melenk & Babuska 1996; Babuska & Melenk 1997). XFEM (Belytschko-Black
1999; Moes-Dolbow-Belytschko 1999) is the local-enrichment specialization
of PU/GFEM applied to crack-growth problems.

The 2021 Jafari et al. paper is, per its own abstract, the **first
documented XFEM implementation inside COMSOL Multiphysics**. It is
notable because COMSOL does NOT expose user-element / UEL hooks the way
ABAQUS does; the enrichment had to be smuggled in through:

  1. several stacked "Solid Mechanics" modules (one per enrichment
     function), each contributing its own DOFs and weak form,
  2. internal-variable surgery via COMSOL's "Equation View" panel
     to override the default strain / stress definitions,
  3. external MATLAB scripts (via LiveLink) for the level-set
     update, geometric search of enriched elements, and crack-tip
     advance after each load step.

## Why XFEM beats remeshing

For a growing / arbitrary crack:

| Approach           | Mesh updates needed | Data transfer | Convergence |
|--------------------|---------------------|---------------|-------------|
| Standard FEM       | yes (after each step) | required      | O(h^p)      |
| Adaptive remeshing | yes (smarter)        | required      | O(h^p)      |
| XFEM               | NO                   | not required  | O(h^p) with proper enrichment |

The mesh stays fixed; the discontinuity moves through the mesh as
enrichment DOFs activate/deactivate on a per-element basis.

## Relation to other enrichment families

| Family              | Year | Scope of enrichment              |
|---------------------|------|----------------------------------|
| PUFEM (Melenk-Babuska)| 1996 | universal, global enrichment    |
| GFEM (Strouboulis et al.)| 2000 | universal, global enrichment |
| XFEM (Belytschko-Black, Moes-Dolbow-Belytschko)| 1999 | LOCAL, only near features |
| ngsxfem (cut-FEM)   | 201x | NGSolve add-on, cut-element view |

XFEM is currently a REFERENCE LAYER for the Sugahara lab: it is not yet
used in production radia-mcp pipelines, but it is the canonical answer
when a future user asks "how do I model a magnetic core with a propagating
fracture in COMSOL?" or "can I avoid remeshing my IH workpiece if it cracks
during heating?".

## Solid-mechanics XFEM vs Electromagnetic XFEM (Hiruma 2023)

THIS module covers **solid-mechanics XFEM** (Jafari et al. 2021, in
COMSOL Multiphysics):  Heaviside + asymptotic crack-tip enrichment,
level-set crack tracking, mixed-mode SIF extraction.

A **structurally distinct** electromagnetic XFEM exists, covered by the
sibling tool ``fem_xfem_em_hiruma``:

| Aspect          | Solid-mech XFEM (this) | EM-XFEM (Hiruma 2023)  |
|-----------------|------------------------|------------------------|
| Reference       | Jafari et al. 2021     | Hiruma et al. IEEE TMag 59(5), 2023 |
| Implementation  | COMSOL + LiveLink/MATLAB | NGSolve compound H1 x H1 |
| Enrichment      | Heaviside + 4 asymptotic | exp(-gamma * xi) (one) |
| Goal            | Stress singularities     | Skin-depth resolution  |
| Domain          | Cracked solid            | Conductor with eddy current |
| Mesh treatment  | Cut by crack             | Conforms to boundary   |
| Lab artifact    | examples/cracked_iron/   | examples/hiruma_xfem_comparison/ |

If the user's question is about eddy currents / induction heating /
high-frequency EM (rather than fracture mechanics), route them to
``fem_xfem_em_hiruma`` instead.

See also:
    - fem_elements('xfem')   - quick XFEM section in the elements catalog
    - fem_potential_formulations('catalog')
    - fem_xfem_em_hiruma     - electromagnetic XFEM sibling
    - radia_mcp.electromagnet (no XFEM yet; possible multi-physics coupling)
"""


THEORY = r"""
# XFEM theory (Jafari et al. 2021, COMSOL formulation)

## 1. Governing equations (linear elastostatics)

A cracked body Omega is bounded by Gamma = Gamma_u U Gamma_t and crack
surface Gamma_d. The equation of motion is

    div(sigma) - rho*u'' + rho*b = 0    in Omega

with boundary conditions

    u = u_bar          on Gamma_u   (Dirichlet)
    sigma * n = t_bar  on Gamma_t   (Neumann)
    sigma * n_Gd = t_d on Gamma_d   (crack-face traction; here zero)

and Hooke's law sigma = D : epsilon (isotropic linear elastic, plane
strain in 2D worked examples; tetra / brick in 3D).

## 2. XFEM displacement approximation

The displacement field is decomposed into a standard + enriched parts:

    u(x) = u_cont(x)
         + M_Gd(x) * u_disc(x)             [Heaviside, jumps across crack]
         + sum_{i=1..4} F_i(x) * u_tip_i(x) [asymptotic tip, 1/sqrt(r)]

where

    M_Gd(x) = H_Gd(phi(x)) - H_Gd(phi(x_I))   (SHIFTED Heaviside, vanishes
                                               at the enriched nodes)
    H_Gd(phi(x)) = +1 if phi(x) >= 0
                 = -1 if phi(x) <  0

and phi(x) is the signed distance from x to the crack interface. The
shift (Liu-Borja 2008) is what makes the enriched DOFs vanish away from
the support of the enrichment, which prevents stiffness-matrix blending
errors at non-enriched neighbours.

The four asymptotic crack-tip functions F_i are the classical Westergaard
near-tip stress fields:

    F_1 = sqrt(r) * sin(theta/2)
    F_2 = sqrt(r) * cos(theta/2)
    F_3 = sqrt(r) * sin(theta/2) * sin(theta)
    F_4 = sqrt(r) * cos(theta/2) * sin(theta)

with (r, theta) the polar coordinates at the crack tip. Only F_1 has the
across-crack jump that captures Mode-I opening; F_2 captures Mode-II.

## 3. Strain field

    epsilon = grad_s u
            = grad_s u_cont
            + H_Gd(phi(x)) * grad_s u_disc
            + delta_Gd * (u_disc kron n_Gd)_s              [Dirac on crack]
            + grad_s ( sum_i F_i(x) u_tip_i )

The Dirac term enforces the discontinuity weakly across Gamma_d. In the
present implementation cracks are traction-free, so the Dirac term drops
out of the weak form (its conjugate boundary integral vanishes).

## 4. Weak form (three coupled blocks)

Test functions are decomposed symmetrically:

    eta = eta_cont + M_Gd(x) * eta_disc + sum_i F_i(x) eta_tip_i

The standard part:

    integral_Omega ( grad_s eta_cont : sigma ) dOmega
      = integral_Omega ( eta_cont . rho b ) dOmega
      - integral_Omega ( eta_cont . rho u'' ) dOmega
      + integral_Gamma_t ( eta_cont . t_bar ) dGamma

The discontinuous (Heaviside) part:

    integral_Omega ( H_Gd(phi(x)) grad_s eta_disc : sigma ) dOmega
      + (boundary integral on Gamma_d, = 0 for traction-free crack)
      = right-hand sides scaled by H_Gd(phi(x))

The tip-enrichment part:

    integral_Omega ( grad_s ( sum_i F_i(x) eta_tip_i ) : sigma ) dOmega
      = right-hand sides scaled by F_i(x)

The support of the second integral is Omega_h (the band of width epsilon
around the crack along which the Heaviside is non-trivial). The support
of the third is Omega_tip (a few elements around the tip).

## 5. Fracture criteria

### Interaction integral for SIFs

To extract K_I, K_II from a numerical solution, the **interaction
integral** is used. Take any closed contour Gamma_J that encloses the
crack tip in the local crack frame x'-y'. The J-integral is

    J = integral_Gamma_J [ W * n_x' - (sigma . d/dx' u) . n ] dGamma

For an actual field "(1)" and an auxiliary near-tip field "(2)" with
pure-K loading:

    I^{1+2} = integral_Gamma_J [ W^{(1,2)} * n_x' - (sigma^{(1)} d/dx' u^{(2)}
              + sigma^{(2)} d/dx' u^{(1)}) ] dGamma

The relation

    I^{1+2} = 2/E' ( K_I^{(1)} K_I^{(2)} + K_II^{(1)} K_II^{(2)} )

with E' = E for plane stress, E' = E/(1-nu^2) for plane strain, gives
K_I^{(1)} by setting (K_I^{(2)} = 1, K_II^{(2)} = 0), and K_II^{(1)} by
setting (K_I^{(2)} = 0, K_II^{(2)} = 1).

In COMSOL this is evaluated using the built-in path-integral operator
`circint(...)` along a circle around the crack tip, centered by the
`at2(...)` coordinate-system operator. A `diskint` (domain form) is the
alternative.

### Maximum hoop stress criterion (Erdogan-Sih)

The crack propagation angle theta_c (measured from the current crack
tangent) is given by

    theta_c = arccos( ( 3 K_II^2 + sqrt(K_I^4 + 8 K_I^2 K_II^2) )
                      / ( K_I^2 + 9 K_II^2 ) )

with the sign of theta_c taken opposite to the sign of K_II. The
equivalent SIF used to test whether propagation occurs is

    K_eq = K_I * cos(theta_c/2)^3
         - 1.5 * K_II * cos(theta_c/2) * sin(theta_c)

Crack propagation is triggered when K_eq > K_IC (material toughness).
A predefined increment dL is added in the global direction (alpha_local +
alpha_global) and the next quasi-static step is solved on the updated
geometry.

The arccos / sign logic is implemented verbatim in
``S:/COMSOL/2022_05_12_XFEM/.../update_crack.m`` (see LAB_FILES).
"""


COMSOL_WORKFLOW = r"""
# COMSOL-specific implementation (Jafari et al. 2021)

COMSOL Multiphysics does NOT have an UEL / user-subroutine entry point
of the same shape as ABAQUS. The 2021 paper proposes a workaround that
exploits three COMSOL features:

  1. multiple Solid Mechanics modules can be stacked, each contributing
     its own DOFs and bilinear-form integral to the global stiffness;
  2. "Equation View" exposes the internal definitions (strain, stress,
     gradient, ...) of every module so they can be overridden in place;
  3. LiveLink for MATLAB executes external .m files as expressions
     within COMSOL, and they may be called during the geometry
     definition, study sequence, and post-processing stages.

## Six Solid Mechanics modules

The standard XFEM split

    u = u_cont + M_Gd * u_disc + sum_{i=1..4} F_i * u_tip_i

is mapped onto SIX Solid Mechanics module instances:

| Module | Symbol     | DOFs                        | Strain definition           |
|--------|------------|-----------------------------|-----------------------------|
| 1      | SM_std     | u_cont (continuous part)    | epsilon = grad_s u_cont     |
| 2      | SM_enr     | u_disc (Heaviside)          | epsilon = H_Gd(phi(x)) grad_s u_disc |
| 3..6   | SM_tip1..4 | u_tip_1..4 (one per F_i)    | epsilon = grad_s (F_i(x) u_tip_i)  |

All six share the same mesh and the same geometry. The stress field is
modified in EVERY module to use the TOTAL strain:

    sigma_total = D : (epsilon_cont + epsilon_disc + epsilon_tip)

so the modules see a consistent stress field even though each contributes
its own strain.

The choice of FOUR tip modules (one per Westergaard function) is what
makes the COMSOL implementation more expensive than a single user-element
in ABAQUS UEL, but it is what permits the proposal to work WITHOUT kernel
modifications to COMSOL.

## Enriched-zone control via Prescribed Displacement

COMSOL does not expose per-node DOF activation. To restrict u_disc and
u_tip_i to their enrichment supports Omega_h and Omega_tip, the paper
exploits the "Prescribed Displacement" domain constraint with a
field-variable indicator psi(x):

    u_disc(x) = u_disc(x)   where psi(x) = 1
              = 0           otherwise

psi is computed in MATLAB from the geometry (1 at enriched nodes, 0
elsewhere) and interpolated to any Gauss point via
``scatteredInterpolant`` in ``interpol.m``. Constraining the unused DOFs
to zero kills both the wasted compute AND the singular block in the
stiffness matrix that would otherwise sit on Omega \ Omega_h.

## MATLAB helper scripts

| Script             | Purpose                                              |
|--------------------|------------------------------------------------------|
| ``preprocess.m``   | Read exported mesh ``mesh.mphtxt`` -> classify which |
|                    | elements / nodes are enriched (geometric search of   |
|                    | crack vs element edges); persist as ``*.mat``        |
| ``phi.m``          | Return the (shifted) Heaviside value at arbitrary    |
|                    | (x, y); supports tip exclusion zone of width 2*epss  |
| ``interpol.m``     | scatteredInterpolant of the enriched-node indicator  |
| ``read_crack.m``   | Return current crack-tip coordinate (complex x+iy)   |
| ``last_angle.m``   | Return angle alpha of the last crack segment         |
| ``update_crack.m`` | Compute theta_c (max hoop stress), K_eq, and if      |
|                    | K_eq > K_IC append a new tip and re-run preprocess   |
| ``p_poly_dist.m``  | Geometric helper: point-to-polyline distance         |

These are called from COMSOL via "Global Variable Probe" + "Variable"
nodes; LiveLink for MATLAB resolves them as expressions.

The path is recovered from ``fileparts(which('<script>.m'))``, so the
entire bundle must sit in the same folder; the mesh export
``mesh.mphtxt`` and the persisted ``*.mat`` files live alongside.

## Numerical integration (CRITICAL)

Standard FEM integrates polynomials with low-order Gauss rules. XFEM
elements crossed by a crack contain a Heaviside discontinuity and need a
much higher integration order, or sub-cell splitting. The paper opts for
RAISING the Gauss order, because COMSOL does NOT expose sub-cell
triangulation:

    SM_enr  : 35-point Gauss quadrature per element
    SM_tip  : 40-point Gauss quadrature per element

These are NGSolve-equivalent "high-order quadrature" choices. The
trade-off is computational cost; the alternative (sub-triangulation) is
not accessible without going into the COMSOL C++ kernel.

A small-support filter (parameter ``delta`` in ``preprocess.m``,
nominally 0.01..0.002) excludes elements where the crack cuts an edge
within ``delta * h`` of a node, to avoid an ill-conditioned support.

## Quasi-static sweep and "fictitious" sweep parameter

For crack propagation, the load is applied incrementally via the
"Auxiliary Sweep" feature ("Stationary" study). COMSOL caches the
geometry from the first sweep value, so geometry-dependent variables
(displacement constraint, modified strain) are forced to be
re-evaluated by adding a fictitious  + lambda * sp  term, with sp the
sweep parameter and lambda ~ 1e-12 (small enough to not perturb the
solution, large enough to invalidate the cache).

## Step-by-step algorithm (Algorithm 1, p.14)

    1. Global Definitions
         - constants (material, load)
         - register MATLAB functions (phi, interpol, read_crack, ...)
    2. Geometry (2D or 3D)
    3. Local Variables
         - interaction integral expressions (W^{1,2}, sigma derivatives)
         - Global Variable Probe calls for crack update / angle / tip
    4. Physics: six Solid Mechanics modules
         (a) SM_std    : material + shape function; modify stress (Eq. 18)
         (b) SM_enr    : material + shape function;
                         modify strain (Eq. 16); modify stress (Eq. 18);
                         apply psi field as Domain Prescribed Displacement
         (c) SM_tip1-4 : material + shape function;
                         modify strain (Eq. 17); modify stress (Eq. 18)
    5. BC + ICs
    6. Mesh
    7. Study: Stationary + Auxiliary Sweep (parameter sp = load fraction)
    8. Post-process: SIF circle integrals, crack trajectory, von Mises
"""


WORKED_EXAMPLES = r"""
# Worked examples in the paper (validated against analytic / experimental)

All examples use linear elastic isotropic material; the 2D cases assume
plane strain. Bi-linear quadrilaterals for 2D; tetra + brick mixture for
3D. Default material: E = 200 GPa, nu = 0.3 (unless stated otherwise).

## 4.1 Center crack in an infinite domain (Section 4.1)

  - Square plate w = 5 m, inclined center crack 2a = 0.2 m -> w/a = 50
    (infinite-domain proxy).
  - Applied: uniaxial 1 MPa traction at top edge, fixed bottom.
  - Analytical SIF: K_I = sigma sqrt(pi a) cos^2(beta),
                    K_II = sigma sqrt(pi a) sin(beta) cos(beta).
  - Mesh refinement at the tip with a/s in {2, 5, 6.7, 9.1, 12.5}.
  - With tip enrichment, error in K_I drops below 0.1% at a/s = 9.1.
  - Energy-norm convergence rate measured at m ~= 0.96 (R^2 = 0.99),
    the theoretically optimal rate for XFEM with shifted Heaviside.

## 4.2 Square plate with 17 randomly distributed cracks

  - 1 m x 1 m, crack length 0.2 m each.
  - 1 MPa tensile traction on top, fixed bottom.
  - XFEM mesh: 14,641 quad elements. FEM reference: 15,000 quads with
    each crack modelled as a finite-width body (0.01 m wide cylinder).
  - Crack opening displacement (COD) agrees to <2% with FEM along an
    arbitrary highlighted crack. Stress field sigma_yy matches
    qualitatively even where multiple tip fields interact.

## 4.3.1 Mixed-mode crack growth in a plate with a hole

  - Geometry from Giner et al. (ABAQUS XFEM benchmark): rectangular plate
    120 x 65 mm with one offset hole; initial crack a_0 = 10 mm; growth
    increment dL = 3 mm.
  - Aluminum alloy: E = 71.7 GPa, nu = 0.33.
  - Mesh: 7,601 quad elements, average size 0.67 mm.
  - The crack curves around the hole; COMSOL trajectory matches the
    experimental fracture path of Giner et al. very closely.
  - **This is the example shipped as the worked .mph + MATLAB bundle**
    in S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with
    hole/. See LAB_FILES["mixed_mode_mph"] and surrounding *.m files.

## 4.3.2 Plate with two holes and two propagating cracks

  - Bouchard et al. (2003) original benchmark.
  - 20 x 10 mm plate, two centered holes 3 mm apart, two pre-cracks of
    length 1 mm; K_IC = 47.4 MPa sqrt(m).
  - Mesh: 12,149 quads (structured).
  - Cracks initially deviate towards the nearer hole then realign with
    their initial direction, matching the adaptive-FEM reference solution
    of Khoei et al. (2008).

## 4.4 Penny-shaped crack in 3D

  - 1 m cube containing a penny-shaped crack of radius 0.25 m at its
    center; 10 MPa tensile traction on top, fixed bottom.
  - Mesh: 50,280 brick elements near crack + 100,843 tet elsewhere
    (average crack-zone size 13 mm).
  - Penny-shape geometry defined ANALYTICALLY in COMSOL (avoiding the
    full level-set MATLAB pipeline for simplicity).
  - FEM reference: same crack modelled as narrow cylindrical void.
  - COD agrees to within 3% with FEM; sigma_zz contours match.
  - Extended further to SIX penny-shaped cracks (radius 0.15 m) parallel
    to the cube faces; results confirm 3D multi-crack capability.
  - **The 3D .mph is shipped** as LAB_FILES["penny_3d_mph"].
"""


CRACK_PROPAGATION = r"""
# Crack propagation: SIF -> direction -> increment

## Pipeline

The propagation step at the end of each quasi-static load increment is:

  1. Solve standard + enriched stress field at current geometry.
  2. Evaluate interaction integrals on a small circle around each crack
     tip -> obtain K_I^{(1)}, K_II^{(1)} for the actual field.
  3. Compute equivalent SIF K_eq using max-hoop-stress (Erdogan-Sih):

        theta_c = sign(-K_II) * acos( ( 3 K_II^2 + sqrt(K_I^4 + 8 K_I^2 K_II^2) )
                                       / ( K_I^2 + 9 K_II^2 ) )

        K_eq    = K_I * cos(theta_c/2)^3
                - 1.5 * K_II * cos(theta_c/2) * sin(theta_c)

  4. If K_eq > K_IC (material toughness):
        a) advance the tip by a predefined dL in the global direction
           alpha_local(theta_c) + alpha_global(previous segment),
        b) append new tip coordinate to the persisted ``crackpt.mat``,
        c) re-run the enrichment classification (preprocess.m logic)
           on the NEW crack geometry,
        d) OR-merge with the existing enrichment indicator psi(x).
     Else: stop, the current load is sub-critical at this tip.

  5. COMSOL re-solves at the next load step on the updated psi field.

## Direction selection criteria (alternatives)

The paper uses the **maximum hoop stress (MTS / Erdogan-Sih)** criterion.
Three families of criteria exist; pick MTS unless the problem demands
otherwise:

| Criterion                      | Driver               | Reference            |
|--------------------------------|----------------------|----------------------|
| Maximum hoop stress (MTS)      | sigma_theta_max      | Erdogan & Sih 1963   |
| Maximum energy release rate    | G_max along theta    | Hussain-Pu-Underwood 1974 |
| Strain energy density (S)      | S_min                | Sih 1974             |
| Maximum principal stress       | sigma_1_max          | (engineering rule)   |

For LEFM under proportional mixed loading, MTS and Gmax usually give
indistinguishable trajectories; MTS is cheaper because it has a
closed-form theta_c.

## Equivalent SIF vs K_IC

The hoop-stress equivalent SIF K_eq blends K_I and K_II:

    K_eq = K_I cos^3(theta_c/2) - (3/2) K_II cos(theta_c/2) sin(theta_c)

A crack advances when K_eq exceeds K_IC. K_IC is a material parameter
(critical mode-I toughness, MPa sqrt(m)). Typical values:
high-strength steel ~50, aluminum 2024-T3 ~26, mild steel ~140,
toughened ceramic ~5.

For purely mode-I loading (K_II = 0), K_eq -> K_I and theta_c -> 0.
"""


LAB_INTEGRATION = r"""
# Sugahara lab integration paths

## Where this knowledge lives

This module sits under ``radia_mcp.fem`` and is served by the
``fem_xfem_comsol`` tool. It is REFERENCE / RESEARCH layer, not a
production pipeline -- there is no ``calc_xfem_*.py`` in
``src/radia/panels/`` at the time of writing.

## When to point a user at this knowledge

  - User asks about fracture mechanics in COMSOL.
  - User asks "can I avoid remeshing my IH workpiece when it cracks?"
  - User asks about Heaviside / asymptotic enrichment in the abstract
    sense (vs cut-FEM, which would point them at ngsxfem).
  - User imports a COMSOL .mph and wants to add level-set crack tracking.

## Related knowledge in radia-mcp

  - ``fem_elements('xfem')`` -- one-paragraph XFEM in the catalog.
  - ``fem_overview('history')`` -- where XFEM sits in the FEM family tree.
  - ``radia_mcp.bem`` -- alternative for stress-singularity problems is
    HDivSurface BEM with analytic crack-tip kernels (different family).
  - ``radia_mcp.electromagnet`` -- a future XFEM+EM coupling for
    piezoelectric or magnetostrictive crack analysis would belong there;
    the paper itself only covers solid mechanics.
  - ``radia_mcp.differential_forms`` -- Heaviside enrichment is a 0-form
    discontinuity; the same machinery in the de Rham complex is what
    enables cut-FEM on HCurl spaces.

## Running the worked example (if the user wants to reproduce)

The mixed-mode-plate-with-hole bundle is self-contained:

```
S:/COMSOL/2022_05_12_XFEM/Mixed-mode crack growth in plate with hole/
    XFEM mixed-mode crack growth in plate with hole.mph
    preprocess.m       <- regenerate crackpt.mat, interpolate.mat, ...
    phi.m              <- called from COMSOL as phi(x, y)
    interpol.m         <- called from COMSOL as interpol(x, y)
    read_crack.m
    last_angle.m
    update_crack.m
    p_poly_dist.m
    mesh.mphtxt        <- pre-exported COMSOL mesh
    crackpt.mat, delta.mat, epss.mat, element.mat,
    node.mat, interpolate.mat  <- regenerated by preprocess.m
```

Per ``READ ME.txt`` in the bundle, the minimal user steps are:

  1. cd into the bundle folder; run ``preprocess.m`` in MATLAB.
  2. Open the .mph in COMSOL with LiveLink for MATLAB enabled.
  3. Compute -> the stationary auxiliary sweep advances the load and
     ``update_crack.m`` propagates the tip whenever K_eq > K_IC.

The 3D penny-shaped example
``S:/COMSOL/2022_05_12_XFEM/XFEM 3D penny-shaped crack.mph`` is
self-contained; the crack geometry is built analytically in COMSOL.

The lab does NOT currently script-drive these .mph files; if a future
panel ever does (e.g. via COMSOL batch CLI), the recommended invocation
is:

```python
import subprocess
from radia_mcp.fem.xfem_comsol_knowledge import LAB_FILES

subprocess.run([
    "comsolbatch",
    "-inputfile", LAB_FILES["mixed_mode_mph"],
    "-outputfile", "C:/temp/xfem_run.mph",
    "-batchlog", "C:/temp/xfem_run.log",
])
```

(COMSOL batch needs a separate license and may require an external
MATLAB to evaluate the LiveLink expressions.)
"""


GLOSSARY = r"""
# XFEM glossary (paper-specific notation)

| Symbol         | Meaning                                                  |
|----------------|----------------------------------------------------------|
| Omega          | Solid body (cracked domain)                              |
| Gamma_u        | Dirichlet boundary                                       |
| Gamma_t        | Neumann (traction) boundary                              |
| Gamma_d        | Crack surface                                            |
| Omega_h        | Heaviside-enrichment support (band of width epsilon)     |
| Omega_tip      | Asymptotic-tip enrichment support (few elements)         |
| u_cont         | Continuous (standard) part of displacement               |
| u_disc         | Discontinuous (Heaviside) enriched DOFs                  |
| u_tip_i        | Asymptotic crack-tip enriched DOFs (i = 1..4)            |
| phi(x)         | Signed distance from x to the crack interface            |
| H_Gd(.)        | Heaviside function (+1 above crack, -1 below)            |
| M_Gd(x)        | SHIFTED Heaviside, M_Gd(x) = H(phi(x)) - H(phi(x_I))     |
| F_i(x)         | Four Westergaard asymptotic crack-tip functions          |
| psi(x)         | Enriched-element indicator field (1 at enriched, else 0) |
| K_I, K_II      | Mode-I and mode-II stress intensity factors              |
| K_IC           | Critical mode-I toughness (material parameter)           |
| K_eq           | Equivalent SIF from max-hoop-stress criterion            |
| theta_c        | Crack propagation angle (relative to current tip frame)  |
| epss           | Heaviside-band half-width (~2 * average element size)    |
| delta          | Small-support filter ratio (typ. 0.01..0.002)            |
| dL             | Predefined crack-tip advance per propagation step        |
| sp             | COMSOL auxiliary sweep parameter (load fraction)         |

## Operator notation

    grad_s u         symmetric gradient (= 0.5 (grad u + grad u^T))
    (.)_s            symmetric part of a tensor
    a kron b         outer product of two vectors -> 2nd order tensor
    delta_Gd         Dirac delta on the crack surface Gamma_d
    E'               = E (plane stress); = E/(1-nu^2) (plane strain)
"""


# ---------------------------------------------------------------------------
# Topic dispatcher
# ---------------------------------------------------------------------------

TOPICS = {
    "overview":          "High-level summary + reference paper context (DEFAULT)",
    "theory":            "Governing equations, PU enrichment, shifted Heaviside, asymptotic tip functions",
    "comsol_workflow":   "Six Solid Mechanics modules + Equation View + MATLAB LiveLink pattern",
    "worked_examples":   "5 paper benchmarks: center crack, random cracks, plate with hole, double hole, 3D penny",
    "crack_propagation": "SIF interaction integral -> theta_c (max hoop stress) -> K_eq vs K_IC -> tip advance",
    "lab_integration":   "Where this fits in radia-mcp; how to run the lab .mph bundle",
    "glossary":          "Notation, symbols, supports, parameters used in the paper",
    "lab_files":         "Absolute paths to the paper PDF, .mph models, MATLAB scripts (LAB_FILES dict)",
    "all":               "Everything (concatenated)",
}


def _lab_files_block() -> str:
    """Render LAB_FILES as a markdown table for display via the MCP tool."""
    lines = [
        "# Lab files (S:/COMSOL/2022_05_12_XFEM/)",
        "",
        "| Key | Path |",
        "|-----|------|",
    ]
    for k, v in LAB_FILES.items():
        lines.append(f"| {k} | {v} |")
    lines.append("")
    lines.append("Use these via:")
    lines.append("    from radia_mcp.fem.xfem_comsol_knowledge import LAB_FILES")
    lines.append("    p = LAB_FILES['mixed_mode_mph']")
    return "\n".join(lines)


def get_xfem_knowledge(topic: str = "overview") -> str:
    """Dispatch XFEM-in-COMSOL knowledge topics.

    Topics:
        overview          - High-level summary (DEFAULT)
        theory            - PU enrichment, shifted Heaviside, asymptotic tip
        comsol_workflow   - Six Solid Mechanics modules + Equation View
        worked_examples   - The five benchmark cases in the paper
        crack_propagation - SIF -> theta_c -> K_eq vs K_IC
        lab_integration   - Where XFEM fits in radia-mcp + reproducer recipe
        glossary          - Notation and symbols
        lab_files         - Absolute paths to PDF / .mph / .m files
        all               - Everything (concatenated)
    """
    t = topic.lower().strip()
    if t in ("overview", "summary", "intro"):
        return OVERVIEW
    if t in ("theory", "math", "formulation"):
        return THEORY
    if t in ("comsol_workflow", "workflow", "comsol", "implementation"):
        return COMSOL_WORKFLOW
    if t in ("worked_examples", "examples", "benchmarks", "verification"):
        return WORKED_EXAMPLES
    if t in ("crack_propagation", "propagation", "sif", "fracture"):
        return CRACK_PROPAGATION
    if t in ("lab_integration", "integration", "lab", "reproduce"):
        return LAB_INTEGRATION
    if t in ("glossary", "notation", "symbols"):
        return GLOSSARY
    if t in ("lab_files", "files", "paths"):
        return _lab_files_block()
    if t == "all":
        return "\n\n".join([
            OVERVIEW,
            THEORY,
            COMSOL_WORKFLOW,
            WORKED_EXAMPLES,
            CRACK_PROPAGATION,
            LAB_INTEGRATION,
            GLOSSARY,
            _lab_files_block(),
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        + ", ".join(TOPICS.keys())
        + "."
    )
