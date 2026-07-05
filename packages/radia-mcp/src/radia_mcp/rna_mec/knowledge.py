"""Reluctance Network Analysis (RNA) / Magnetic Equivalent Circuit (MEC).

Knowledge base derived from canonical RNA / MEC references:

  - Derbas, Williams, Koenig, Pekarek, "A Comparison of Nodal- and
    Mesh-Based Magnetic Equivalent Circuit Models", IEEE Trans. Energy
    Conversion 24(2), 388-395, 2009.
  - Janet, Coulomb, Chillet, Mas, "Magnetic Moment and Reluctance
    Network Mixed Method Applied to Transformer's Modeling", IEEE
    Trans. Magn. 41(5), 1428-1431, 2005.
  - Janet, Coulomb, Chillet, Mas, "Simplified Magnetic Moment Method
    Applied to Current Transformer Modeling", IEEE Trans. Magn. 40(2),
    818-821, 2004.
  - Hane, Nakamura, "Dynamic Hysteresis Modeling for Magnetic Circuit
    Analysis by Incorporating Play Model and Cauer's Equivalent Circuit
    Theory", IEEE Trans. Magn., DOI 10.1109/TMAG.2020.3004355, 2020.
  - Yin, Naidjate, Bracikowski, Pierquin, Trichet, "Topology
    Optimization of Magnetic Actuator based on Reluctance Network
    Modeling and Adjoint Variable Method", Compumag/CEFC, 2023.
  - Lee, Lee, Choi, Park, "Reduced Modeling of Eddy Current-Driven
    Electromechanical System Using Conductor Segmentation and Circuit
    Parameters Extracted by FEA", IEEE Trans. Magn. 41(5), 1448-1451,
    2005. (TEAM Workshop Problem 28)
  - Kameari, Ebrahimi, Sugahara, Shindo, Matsuo, "Cauer Ladder Network
    Representation of Eddy-Current Fields for Model Order Reduction
    Using Finite-Element Method", IEEE Trans. Magn. 54(3), 7201804,
    2018.
  - Le-Duc, Chadebec, Guichon, Meunier, Lembeey, "Coupling between
    PEEC and magnetic moment method", COMPEL 32(1), 383-395, 2013.

All knowledge below is plain ASCII; no Unicode mathematics symbols are
used.  Coordinates in meters, B in Tesla, H in A/m (Radia convention).
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "overview": "RNA / MEC landscape, lab lineage, when to use vs FEM",
    "mec_basics": "Reluctance, permeance, MMF, Hopkinson's law, BH branch",
    "nodal_vs_mesh_analysis": "Derbas 2009: nodal (KCL) vs mesh (KVL), Jacobian conditioning",
    "reluctance_network_construction": "Flux tube discretization, claw-pole example, area selection",
    "lumped_extraction_fea": "Inductance L/M extraction from FEM, conductor segmentation",
    "cauer_ladder_rna": "CLN representation of eddy-current fields, Kameari 2018, Hane play+Cauer",
    "rna_magnetic_coupling": "RNA + magnetic-moment mixed method, Janet transformer, current transformer",
    "electromechanical_coupling": "TEAM-28 levitation, Runge-Kutta ODE, dL/dz force",
    "team28_reduced_model": "TEAM-28 reduced model in depth, 85h FEM -> 1h CLN+segmentation",
    "topology_optimization": "RNA + AVM (adjoint variable method), SIMP, magnetic actuator TO",
    "dynamic_hysteresis": "Play + Cauer dynamic hysteresis MEC (lab specialty)",
    "vs_pec_peec": "RNA vs PEEC vs PEC (partial element equivalent circuit) terminology",
    "all": "All topics concatenated",
}


OVERVIEW = r"""
# Reluctance Network Analysis (RNA) / Magnetic Equivalent Circuit (MEC)

RNA (also called MEC = Magnetic Equivalent Circuit) discretizes a
magnetic structure into a network of lumped magnetic elements --
reluctances (or permeances), MMF sources, and flux sources -- and
solves it like an electric circuit using Kirchhoff laws.

## The magnetic-electric analogy

| Electric quantity      | Magnetic counterpart                       |
|------------------------|--------------------------------------------|
| Voltage V              | MMF F = N I (or H integrated along path)  |
| Current I              | Flux Phi = integral B dA                  |
| Resistance R           | Reluctance R = l / (mu A)                 |
| Conductance G          | Permeance P = mu A / l = 1/R              |
| Ohm's law V = R I      | Hopkinson's law F = R Phi                 |
| Kirchhoff current law  | sum of branch fluxes at a node = 0        |
| Kirchhoff voltage law  | sum of MMF drops around a loop = sum F_src|

The analogy is exact for LINEAR magnetostatics in a tube-of-flux
approximation.  For nonlinear iron (BH curve), the analogy is local:
reluctance R(Phi) = l / (mu(B) A), where B = Phi / A, mu = mu(B), and
the network becomes a nonlinear algebraic system.

## When to use RNA / MEC versus FEM

| Scenario                              | RNA / MEC        | FEM            |
|---------------------------------------|------------------|----------------|
| Quick concept design                  | best             | overkill       |
| Lumped equivalent for SPICE coupling  | natural          | hard           |
| Detailed flux distribution            | rough            | best           |
| Saturation in core                    | OK (lookup BH)   | OK (BH curve)  |
| Hysteresis (MEC + Play / Energy)      | LAB specialty    | OK             |
| Leakage flux (3D)                     | difficult        | natural        |
| Optimization (sensitivity / TO)       | fast             | slow           |
| Population-based design (10^6 calls)  | only choice      | infeasible     |
| Frequency sweep with eddy currents    | Cauer ladder OK  | per-freq FEM   |
| Multi-physics with mechanical motion  | natural (ODE)    | re-mesh per dt |

RNA loses spatial detail; FEM loses speed and lumped-port intuition.
In practice they cohabit: use FEM as the trusted reference (and the
parameter extractor), use RNA as the dynamic / optimization model.

## Three flavours of RNA / MEC

1. **Pure RNA**: human picks the reluctance topology by inspection
   (Roters 1941, Slemon 1953, Ostovic 1986).  Few branches, fast.
2. **RNA from FEA extraction**: lumped parameters R, L, M derived from
   per-segment FEM computations (Lee 2005 TEAM-28).  Bridges FEM
   accuracy with RNA speed.
3. **Mesh-based RNA** (regular grid): mesh the design region with a
   grid of reluctance cells, each cell = 4 reluctances + flux sources.
   Used in topology optimization (Yin 2023).

## Lab core papers (public-safe curated corpus)

Tanaka / Hane lineage (Tohoku University Nakamura lab, collaborating
with Sugahara lab):
  - `1609_マグ研_田中.pdf` -- Tanaka 2016
  - `1702_マグ研_田中.pdf` -- Tanaka 2017
  - `1801_磁気学会_羽根_採録.pdf` -- Hane 2018 (Magnetics Society of Japan)
  - `1803_マグ研_羽根.pdf` -- Hane 2018
  - `1902_マグ研_羽根.pdf` -- Hane 2019
  - `1912_マグ研_羽根_PlayCauer.pdf` -- Hane 2019, Play + Cauer
  - `修士論文_本審査_久田.ppt` -- Hisada (M.Eng. thesis)
  - `博士論文_本審査_吉田.pptx` -- Yoshida (Ph.D. thesis)
  - `AC_DC_Bridgeless_Flyback.pdf` -- application

External (`public-safe curated corpus`):
  - `Comparison_Nodal_Mesh_MEC.pdf` (Derbas 2009) -- the canonical
    nodal-vs-mesh comparison.
  - `Dynamic_Hysteresis_Magnetic_Circuit.pdf` (Hane 2020) -- the
    English IEEE version of the Play + Cauer MEC.
  - `TopOpt_Magnetic_Actuator_RNA.pdf` (Yin 2023) -- RNA + Adjoint
    Variable Method for topology optimization.

Hysteresis-focused (`public-safe curated corpus`):
  - `(1)LSモデル(Physica B, 2000).pdf` -- LS hysteresis model
  - `(2)LSモデル_RNA(Not peer-reviewed,2017).pdf` -- LS embedded in RNA
  - `プライザッハRNA(東南大学,IEEE,2019).pdf` -- Preisach RNA (Southeast
    University, China)

## Cross-references in radia_mcp

- `radia_mcp.radia_ngsolve.hdiv_vim` -- magnetic-material field coupling
- `radia_mcp.magnetic_materials.hysteresis_models.play` -- Play model
- `radia_mcp.magnetic_materials.hysteresis_models.lab_core` -- Energy /
  Play core (Tohoku-Kindai lineage)
- `radia_mcp.mor.systematic.cln` -- Cauer Ladder Network (CLN) MOR
- `radia_mcp.motor.calc_motor_transient` -- time-stepping framework
"""


MEC_BASICS = r"""
# MEC basics: reluctance, permeance, MMF, BH branch

## Hopkinson's law (the magnetic Ohm's law)

For a homogeneous flux tube of length l, uniform cross-section A and
permeability mu:

    R = l / (mu A)       [A-turn / Wb]   (reluctance)
    P = mu A / l = 1/R   [Wb / A-turn]   (permeance)
    F = R Phi            (Hopkinson's law, F in A-turn, Phi in Wb)

R is the magnetic equivalent of electrical resistance; F (MMF) is the
equivalent of voltage; Phi (flux) of current.

## Three branch types

1. **Linear reluctance** (air, linear iron approximation):
   R = const, mu = mu_0 mu_r constant.  In an air gap of length g and
   area A_gap: R_gap = g / (mu_0 A_gap).

2. **Nonlinear (BH curve)** branch (saturating iron):
   R = R(Phi) = l / (mu(B) A) where B = Phi / A and mu(B) is obtained
   from a measured BH or mu(B) curve.  Common implementations use a
   piecewise quadratic spline (Derbas 2009 Fig. 4-5) or a Frohlich /
   Brauer parametrisation.

3. **MMF source** (winding): F_src = N I where N is the number of
   turns and I the winding current.  A permanent magnet contributes a
   constant MMF F_pm = H_c t_mag where H_c is coercivity and t_mag the
   magnet thickness, in series with the magnet's internal permeance
   P_m = mu_0 mu_rmag l_mag / t_mag (Derbas 2009 eqs. 56-57).

## Flux-tube discretization in 2D / 3D

The "flux tube" is the central assumption: pick a region where flux is
roughly axial, parameterize as l (length along flux) and A (cross
sectional area perpendicular to flux).  When the cross-section varies
along the tube (e.g. trapezoidal claw, conical pole tip), integrate:

    R = integral over tube of  dz / (mu(z) A(z))                   (1)

For a linearly-tapered claw section of widths l_1 (entry) and l_2
(exit) and length l_claw (Derbas 2009 eq. 61):

    R_ci = l_claw / (mu_0 mu_rci N (l_1 - l_2)) * ln(l_1 / l_2)    (2)

The factor 1 / (l_1 - l_2) appears because dA/dz = const for the
trapezoidal taper; the logarithm absorbs the variable cross section.

## "Chosen area" hazard (Derbas 2009)

For computing B in a tapered tube one needs A_c to convert Phi to B:

    B = Phi / A_c                                                    (3)

The choice of A_c (entry, exit, midpoint, two-thirds-from-base, ...)
is NOT unique.  Derbas reports that the "two-thirds-from-the-longer-
side" choice gives the best match to FEM for the claw-pole example.
This is a known modeling-art residue: RNA accuracy on tapered geometry
depends on a heuristic area selection, and the heuristic is not
universally valid.

## Branch matrix assembly

Once each branch i has (R_i, F_src_i), the full network is assembled
just like a resistive electric circuit:

    Nodal (KCL):  G u = phi_src   where G_ii = sum of permeances
                                  incident on node i and G_ij = - P_ij
                                  for the branch joining i and j.

    Mesh (KVL):   A_R phi_loop = F_src_loop  where A_R is the mesh
                                  reluctance matrix obtained by KVL
                                  around each independent loop.

For LINEAR systems the two give identical results (up to the linear
transform between node potentials u and loop fluxes phi_loop).  For
NONLINEAR (saturating) systems they have DIFFERENT numerical
behaviour -- see topic `nodal_vs_mesh_analysis`.

## What MEC cannot do (without help)

- **3D leakage with no obvious flux tube** -- e.g. quick stray field
  around the corner of a transformer.  Needs RNA + magnetic-moment model or RNA + FEM.
- **Eddy currents in a conductor of arbitrary cross-section** --
  needs Cauer ladder (Kameari 2018) or per-segment lumped extraction
  from FEM (Lee 2005, TEAM-28).
- **Detailed force distribution on a moving armature** -- needs FEM or
  a fine RNA mesh with Maxwell-stress tensor evaluation.

These limitations motivate the four hybrid approaches documented in
the other topics: nodal-vs-mesh, RNA + FEA extraction, RNA + magnetic-moment model
coupling, and RNA + Cauer ladder.
"""


NODAL_VS_MESH_ANALYSIS = r"""
# Nodal vs mesh MEC analysis (Derbas 2009)

The single most useful operational paper on this subject is:

  Derbas, Williams, Koenig, Pekarek, "A Comparison of Nodal- and
  Mesh-Based Magnetic Equivalent Circuit Models", IEEE Trans. Energy
  Conv. 24(2), 388-395, June 2009.

The finding: **under nonlinear / saturated operation, mesh-based MEC
converges roughly an order of magnitude faster than nodal-based MEC**,
because the mesh-based Jacobian has a much smaller condition number.

## The two formulations

**Nodal-based MEC** (Kirchhoff Current Law).  Unknowns = nodal MMF
potentials u.  Permeances P relate node potentials to branch fluxes:

    A_P u = phi_src                                                  (1)

where A_P is the nodal permeance matrix (analogue of nodal admittance
matrix in electrical circuits), and phi_src injects MMF-derived flux
sources at PM nodes.  In the linear case A_P is constant; in the
nonlinear case A_P = A_P(u) since each P_ci depends on mu_rci, which
depends on B_ci, which depends on Phi_ci = P_ci (u_i - u_j).

**Mesh-based MEC** (Kirchhoff Voltage Law).  Unknowns = loop fluxes
phi_loop.  Reluctances R relate branch fluxes to branch MMF drops:

    A_R phi = F_src                                                  (2)

where A_R is the mesh reluctance matrix (analogue of mesh impedance
matrix).  Again, in linear cases A_R is constant; in nonlinear cases
A_R = A_R(phi) because R_ci depends on mu_rci(B_ci) and the branch
flux is a linear combination of loop fluxes.

## Newton-Raphson solution

Both formulations need a nonlinear solve:

    Nodal:    f(u) = A_P(u) u - phi_src = 0,     J = df/du
    Mesh:     g(phi) = A_R(phi) phi - F_src = 0, J = dg/dphi

Apply Newton-Raphson: u_{n+1} = u_n - J^-1 f(u_n) (eq. 55 Derbas).

The Jacobian has two parts: the constant matrix (A_P or A_R) PLUS a
correction matrix capturing the dependence of permeance / reluctance
on the working variable.

## Why the nodal Jacobian is ill-conditioned

For the nodal formulation, the partial derivative chain rule gives:

    dP_c1 / du_m1 = (dP_c1 / dmu_rc1) * (dmu_rc1 / dB_c1)
                  * (dB_c1 / dPhi_c1) * (dPhi_c1 / du_m1)            (3)

The trailing term dPhi_c1 / du_m1 itself depends on dP_c1 / du_m1
(eq. 24 in Derbas), so a self-consistent algebraic rearrangement is
needed (eq. 25 in Derbas):

    dP_c1/du_m1 = (alpha beta gamma P_c1) / (1 - alpha beta gamma * du)
                                                                     (4)

This formula has a singularity when (alpha beta gamma * du) -> 1,
causing the Jacobian to become ill-conditioned and Newton-Raphson to
either diverge or oscillate.

For the mesh formulation, the analogous chain rule (eq. 31 Derbas)
ends in dPhi_c12 / dphi_m = -1, a CONSTANT.  No self-consistent
algebraic rearrangement is needed.  The Jacobian is well-behaved.

## Empirical Newton-Raphson convergence (Derbas Tables I-IV)

Test problem: 2D claw-pole alternator with 3 claw sections,
permanent magnet coercivity H_c swept from 1e5 to 1.8e6 A/m.

| H_c              | Nodal (analytic J) | Nodal (J = A_P) | Mesh    |
|------------------|--------------------|-----------------|---------|
| 1e5 (linear)     | converges          | converges       | converges |
| 8e5 (saturated)  | DIVERGES           | ~50 iter        | ~5 iter |
| 16e5             | DIVERGES           | ~100 iter       | ~7 iter |
| 18e5             | DIVERGES           | ~150 iter       | ~9 iter |

The condition number of the mesh Jacobian stays O(10^2 - 10^3); the
nodal-analytic Jacobian condition number grows to O(10^8 - 10^9) as
the iron saturates further.

## Why scaling does not save the nodal formulation

Derbas considered rescaling the unknown vector u.  The ill-
conditioning is caused by the air gap permeance P_ag being roughly
1000x smaller than iron permeances; P_ag appears in many off-diagonal
entries.  Rescaling one variable that multiplies P_ag forces other
permeances to be rescaled too, and the resulting matrix is again
ill-conditioned.  The structural problem cannot be fixed by scaling.

## Practical recommendation

For nonlinear MEC, **default to mesh-based formulation**:

1. Pick independent loops (often topology trees / co-trees).
2. Loop variables = loop fluxes.
3. Newton-Raphson with analytic Jacobian: a few iterations to converge.

The nodal formulation is fine for LINEAR MEC, or for cases where the
permeances do not depend strongly on the working variable (small-
signal AC at fixed bias).  But the moment iron saturation is
significant, switch to mesh formulation.

## Multiply-connected geometry: cohomology cuts

A NODE-based formulation needs no cuts; the difficulty is purely
numerical (above).  A MESH-based formulation needs a tree-cotree
decomposition (or cohomology cuts) to identify the independent loops
in multiply-connected geometry.  For simple ladder-like cores (claw-
pole, EI-core, three-leg transformer) the cuts are obvious; for
complex 3D structures they require algorithmic cohomology
(see Dlotko / Specogna for general algorithms).
"""


RELUCTANCE_NETWORK_CONSTRUCTION = r"""
# Reluctance network construction: from geometry to A_R

## The four-step process

1. **Identify flux paths**.  By inspection or by an a-priori coarse
   FEM, mark out the dominant flux tubes (legs, yokes, air gaps,
   leakage paths).
2. **Section each tube**.  Subdivide along the flux direction into N
   sections such that the cross-section can be treated as piecewise
   uniform.  For a claw-pole, N = 30 sections per claw is enough
   (Derbas 2009 Fig. 6).
3. **Compute branch parameters**.  Each section i becomes a reluctance
   R_i = l_i / (mu(B_i) A_i) plus optional MMF source F_src_i.
   Air gaps use mu = mu_0; iron sections use mu = mu(B_i) interpolated
   from the BH curve.
4. **Assemble**.  Apply KCL at each node (nodal) OR KVL around each
   loop (mesh) to build A_P or A_R.

## Worked example: claw-pole alternator (Derbas 2009 Appendix)

Geometry: two trapezoidal claws facing each other across an air gap,
with an interior stator core in between.  Each claw is divided into N
sections.  The corresponding air-gap is sectioned into N sub-gaps,
each connected only to its facing claw section.

For section i of claw with widths l_1 (entry, base) and l_2 (exit,
tip), depth z = 1, and N sections along the claw of total length
l_claw, the section length is dz = l_claw / N and the width varies
linearly along z.  Section reluctance:

    R_ci = l_claw / (mu_0 mu_rci N (l_1 - l_2)) ln(l_1 / l_2)        (1)

(Derbas eq. 61.)  The N divides the section reluctance into N pieces
since the tube is subdivided.  The logarithm comes from the analytic
integration of dR = dz / (mu_0 mu_rci A(z)) over the trapezoid.

For the air gap (each sub-section, area l_claw/N x depth, length
l_gap = const along the radial direction):

    R_ag = N l_gap / (mu_0 l_claw)                                   (2)

For the stator (uniform, constant cross section):

    R_s = l_stat / (mu_0 mu_s l_claw)                                (3)

For the permanent magnet (radial bias, coercivity H_c, thickness
t_mag, area l_mag x depth):

    R_m = t_mag / (mu_0 mu_rmag l_mag)                               (4)
    F_pm = H_c t_mag                                                 (5)

## BH curve handling

Derbas builds a piecewise interpolant:
  - Linear segment for small B (below the knee)
  - Quadratic spline through the knee region
  - Constant mu = 1 (saturation, mu_r -> 1) for large B

The mu(B) curve (not mu(H)!) is the working object because the MEC
unknown is Phi (hence B) -- not H.  At each Newton iteration, you
need both mu and dmu/dB (for the Jacobian).  Numerical differentiation
of a noisy mu(B) measurement is the leading source of nonconvergence;
hand-fit a smooth analytic form (Frohlich, Brauer) wherever possible.

## Grid-based RNA (Yin 2023 topology optimization)

For TO, the geometry is unknown and the design region is meshed with a
regular grid of (Nx x Ny) cells.  Each cell has 4 reluctances around
its perimeter and 2 (or 4) flux sources driven by the local current
density:

    R_l = R_r = (1 / (mu_0 mu_r)) (1/z) (dx/2) / dy                  (6)
    R_u = R_d = (1 / (mu_0 mu_r)) (1/z) (dy/2) / dx                  (7)
    phi_l = phi_r = N I / R_l                                        (8)

where dx, dy are the grid step sizes, z is the device depth, and N is
chosen to respect Ampere's law in each loop.  The whole grid network
is assembled as a single matrix [Lambda] eps = phi (eq. 4 Yin), with
[Lambda] containing all the reluctances and eps the magnetic-potential
vector.  This formulation is the magnetic analog of a resistor-grid
electrical network and is the foundation of RNA-based topology
optimization (next topic).

## Connecting RNA to circuit simulators

Each branch (R, F_src) is a SPICE-compatible element:
  - R becomes a passive element in the SPICE netlist
  - F_src = N I becomes a controlled source whose value tracks the
    winding current
  - flux Phi flowing through a branch maps to "current" in the SPICE
    representation; voltage in SPICE maps to MMF.

This gives an immediate Spice-friendly "magnetic-circuit" subnetlist
that can be embedded in a power-electronics simulator (LTspice,
ngspice) to study transformer + drive interaction at low cost.  This
is the classic application of RNA in power electronics design.

## Validating an RNA: comparison to a one-shot FEM

Standard practice: build a coarse FEM model, compute one operating
point (peak DC bias), compare Phi_branch from FEM (by integrating B.n
over the surface of each tube) to Phi_branch from RNA.  Branches that
disagree by more than 5-10% point to:
  - missing leakage paths (add a parallel reluctance)
  - wrong "chosen area" A_c in tapered sections (re-pick A_c)
  - missing fringing (extend air-gap permeance with Schwarz-Christoffel
    correction)

Iterate the RNA topology until agreement is within tolerance, then
trust it for sweep / optimization studies.
"""


LUMPED_EXTRACTION_FEA = r"""
# Lumped parameter extraction from FEA

The reluctance / inductance / mutual-inductance values that populate
an RNA can be PRECOMPUTED by a finite-element solver and then frozen
as constants in the network solve.  This is the bridge between
high-accuracy electromagnetic field analysis (FEM) and high-speed
RNA / circuit simulation.

Two canonical sources:

  - Lee, Lee, Choi, Park, "Reduced Modeling of Eddy Current-Driven
    Electromechanical System Using Conductor Segmentation and Circuit
    Parameters Extracted by FEA", IEEE Trans. Magn. 41(5), 1448-1451,
    2005.  (TEAM Workshop Problem 28.)
  - Kameari et al. 2018 -- CLN extraction (covered in `cauer_ladder_rna`).

## The Lee 2005 procedure (TEAM-28)

TEAM Workshop Problem 28 is a coaxial pair of coils underneath a thin
aluminum plate.  When AC is applied, eddy currents in the plate
oppose the field and the plate is levitated.  The plate's motion
follows a nonlinear ODE coupling position z, velocity v, and
multiple eddy-current loops.

Outline:

1. **Segment the conductor** in the direction of expected eddy current
   loops.  Plate divided into N_seg = 1, 2, or 5 concentric rings.
   Each segment must be thinner than the skin depth so the per-segment
   current can be treated as uniform.

2. **Sweep the moving part over its range**.  The plate's height z is
   varied across the expected motion range; at each height a magneto-
   static FEM is solved.

3. **At each height, extract**:
   - Self-inductance L_kk of each segment k
   - Mutual inductance M_kj between every pair of segments
   - Mutual inductance M_kfixed between each segment and the fixed coil
   - All these are FUNCTIONS of z (since geometry changes as plate
     moves) and stored as discrete tables.

4. **Interpolate**.  Lee uses cubic-spline interpolation of L(z),
   M(z) so that derivatives L'(z) = dL/dz are smooth and continuous.
   The derivatives are essential for the force computation
   F = 0.5 * sum_{ij} I_i I_j dM_ij/dz (eq. 5 Lee).

5. **Assemble state equations**.  State = (z, v, {I_k}).  Equations of
   motion + per-segment Kirchhoff voltage:

       m dv/dt = F_em - m g                                          (1)
       (R_k + R_k_ext) I_k + d(lambda_k) / dt = V_k                  (2)

   where lambda_k = sum_j L_kj I_j is the flux linkage of segment k,
   d(lambda_k)/dt has both an inductance-change term and a current-
   change term:

       d(lambda_k)/dt = sum_j L_kj dI_j/dt
                      + v * sum_j (dL_kj/dz) I_j                     (3)

   The second term is the motional EMF.

6. **Time-step**.  Use a stiff ODE solver (Runge-Kutta 4 with adaptive
   step works; Lee uses RK4 explicit).  At each step, the state
   variables are updated using the pre-computed L(z), M(z), and their
   derivatives -- NO FEM IS CALLED INSIDE THE TIME LOOP.

## Speedup vs full transient FEM

For TEAM-28 (Lee 2005 Table II):

| Method                          | CPU time | Plate height error |
|---------------------------------|---------:|-------------------:|
| Time-stepping FEM (reference)   |     85 h |     reference      |
| Reduced model, 1 segment        |    <1 h  |     14.1%          |
| Reduced model, 2 segments       |    <1 h  |      5.5%          |
| Reduced model, 5 segments       |    <1 h  |    <5.5% (best)    |
| Experiment                      |        - |     reference      |

Two orders of magnitude faster (85 h -> < 1 h), with 5% accuracy at 5
segments.  The reduced model also allows the operator to change the
input voltage, frequency, or initial position WITHOUT re-running FEM:
the L(z), M(z) tables stay valid for any operating condition.  This
is critical for control design and Monte-Carlo studies, where the
full-FEM cost would be prohibitive.

## Conductor segmentation guidance

- Segments must be thinner than the skin depth delta = sqrt(2 / (mu w sigma))
  at the highest frequency of interest, otherwise the assumption of
  uniform current per segment fails.
- Segments should follow the topology of the expected eddy-current
  paths (concentric rings for an axisymmetric plate, longitudinal
  strips for a moving rail, etc.).
- More segments improve accuracy at high frequency but also expand
  the L matrix as N_seg^2.  Practical sweet spot: 3-10 segments.

## Range of validity

The Lee approach assumes:
  - Linear material (no saturation).  If the iron path saturates, the
    L matrix becomes state-dependent and the per-(z) table must be
    extended to per-(z, I_bias).
  - Slow motion compared to the electromagnetic timescale.  At very
    high speeds the FEM extraction at static z is no longer accurate
    (eddy-current drag becomes a strong function of dz/dt).
  - 1 / time-step >> 1 / electromagnetic time constant.  The ODE
    solver must resolve the smallest L/R time constant in the network.

For levitation / motor / actuator transients in the 100 Hz - 10 kHz
range, the assumptions are usually safe.

## Generalisation: parameter extraction for any RNA branch

The procedure is general: any RNA branch's parameter can be extracted
by:
  1. Running an FEM solve with that branch "alone" (other branches
     deactivated / replaced by perfect short / perfect open).
  2. Computing the desired lumped quantity by integration of the
     resulting field (flux, energy, force, ...).
  3. Repeating across the design parameter (current, position,
     saturation level) and tabulating.

The extraction itself is parallel-embarrassing: each operating point
is independent.  This is the bedrock of FEA-extracted RNA libraries
for power electronics and electric machine design.
"""


CAUER_LADDER_RNA = r"""
# Cauer Ladder Network (CLN) representation of eddy-current fields

The Cauer Ladder Network (CLN) is a continued-fraction expansion of a
linear eddy-current impedance Z(s) into a ladder of alternating
resistors and inductors:

    Z(s) = R_0 + 1 / ( L_1 s + 1 / ( R_2 + 1 / ( L_3 s + ... )))     (1)

This represents the input impedance of a 1-port eddy-current device
as a low-order RL ladder, exactly to within truncation error.  CLN
turns out to be a DRAMATIC model-order reduction tool for both 1D
laminated sheets and full 3D FEM eddy-current problems.

The seminal extension to 3D FEM is:

  Kameari, Ebrahimi, Sugahara, Shindo, Matsuo, "Cauer Ladder Network
  Representation of Eddy-Current Fields for Model Order Reduction
  Using Finite-Element Method", IEEE Trans. Magn. 54(3), 7201804,
  2018.

## CLN as MOR for full 3D FEM

The Maxwell eddy-current equations (rotH = sigma E, rotE = -d(muH)/dt)
are decomposed into a sequence of static E-mode and static H-mode
fields, with weights e_{2n}(t) and h_{2n+1}(t) governed by Kirchhoff
laws of the CLN:

    E = sum_n e_{2n}(t) E_{2n}                                       (2)
    H = sum_n h_{2n+1}(t) H_{2n+1}                                   (3)

The modes are extracted recursively:

    Step 0:  Solve static E_0 (1 V applied), compute R_0 = 1/(sigma E_0)^2
    Step 1:  Solve magnetic mode H_{2n-1} from  rot(H_{2n-1}) = R_{2n-2} sigma E_{2n-2}
             Compute L_{2n-1} = mu H_{2n-1}^2
    Step 2:  Solve electric mode E_{2n} from  rot(E_{2n}) = -(1/L_{2n-1}) mu H_{2n-1}
             Compute R_{2n} = 1 / (sigma E_{2n})^2
    Step 3:  Solve CLN in frequency or time
    Step 4:  Reconstruct E, H by superposition

This requires only N (a few) STATIC FEM solves and gives the full
frequency response over a wide band, plus the full time response under
arbitrary excitation, with no per-frequency solve and no per-time-
step solve.

## Performance (Kameari 2018 Table II)

3D nonlinear DC-DC converter inductor, 104k elements, 109k nodes:

| Method                            | CPU time          | Speedup   |
|-----------------------------------|-------------------|-----------|
| Full transient FEM (800 steps)    | ~hours            | 1x        |
| CLN, 1 stage                      | ~tens of seconds  | ~100x     |
| CLN, 5 stages (high freq)         | ~minute           | ~50x      |

Even 1-stage CLN gives accurate steady-state current; 2-stage CLN is
indistinguishable from FEM up to 1 MHz.

## CLN + Play hysteresis for PWM iron loss (Hane 2020)

Hane, Nakamura, "Dynamic Hysteresis Modeling for Magnetic Circuit
Analysis by Incorporating Play Model and Cauer's Equivalent Circuit
Theory", IEEE Trans. Magn., DOI 10.1109/TMAG.2020.3004355, 2020.

The PWM-driven minor-loop iron loss problem combines TWO effects:
  1. Hysteresis (modeled by Play model)
  2. Eddy-current skin effect inside the laminated sheet

The Play model alone is fine at the fundamental but misses the
HF carrier minor loops (which see a shrunken eddy-current path
because the skin depth is now < lamination thickness).  Hane's
solution:

    Replace the single eddy-current inductance L_1 in the previous
    Play-based MEC with a Cauer-I ladder (R_0, L_1, R_2, L_3, ...)
    derived from the sheet's 1D cell problem.

Cauer parameters (Hane eq. 2-3):

    L_2n+1 = mu * geometric_factor                                   (1)
    R_2n   = 4 / (sigma d^2)                                         (2)

where d is the lamination thickness, sigma the conductivity, and mu
the differential permeability evaluated on the centerline of the DC
hysteresis loop at B_m = 1 T (Hane uses the dc-1T BH operating point
as the linearization point for the eddy-current ladder; this is a key
empirical recipe in the paper).

The ladder is truncated at 2 stages for the practical PWM-driven
silicon-steel ring core; this captures the skin effect at 1-2 kHz
carrier with adequate accuracy.  Hane shows experimental validation
on grain-oriented and non-oriented Si-Fe rings.

## Combined Play + Cauer + RNA recipe (Hane 1912 lab paper)

The Tohoku-Kindai "1912" paper extends Hane 2020 by embedding the
Play + Cauer MEC into a multi-branch RNA for a permanent-magnet
motor / variable inductor.  Per-branch structure:

    Branch i = [ Play hysterons + Cauer ladder ] in series

Each branch's instantaneous BH state evolves via:
  - Play model O(K) forward evaluation (B -> H), K play operators
  - Cauer ladder feeding the skin-effect-corrected eddy current
  - Anomalous eddy-current loss as a separate dependent source

In the time loop:
  1. Update winding currents (from external SPICE circuit)
  2. Solve nonlinear KCL with current B-H tangents
  3. Update Play state on every branch
  4. Update Cauer state on every branch
  5. Compute total flux + per-branch loss breakdown (hysteresis +
     classical eddy + anomalous eddy)

This is the **complete dynamic hysteresis MEC** of the Sugahara-
Nakamura-Hane lineage, and is the most physically detailed RNA
available for PWM iron-loss prediction.

## Cross-references

- `radia_mcp.mor.systematic.cln` -- the general CLN MOR knowledge for
  arbitrary 3D FEM problems (Kameari 2018 and successors).
- `radia_mcp.magnetic_materials.hysteresis_models.play` -- Play model
  forward / inverse evaluation.
- public-safe curated corpus -- 30+ Cauer / CLN
  papers (Sato/Igarashi, Matsuo, Shindo, Eskandari, etc.).
"""


RNA_MAGNETIC_COUPLING = r"""
# RNA + magnetic-moment mixed method (Janet 2004, 2005)

Pure RNA struggles with 3D LEAKAGE flux outside the dominant flux
tubes.  Pure magnetic-moment method handles 3D leakage but needs
many elements to resolve a long, thin closed magnetic circuit (slow
mesh refinement, slow MoM matrix assembly).

The Janet (Schneider Electric / LEG Grenoble) work shows how to
**combine** them for current-transformer-class problems: use a few-
element magnetic-moment model to capture 3D leakage + saturation, then ADJUST its
coupling coefficients using a fast RNA computation that captures the
linear-regime behaviour accurately.

Two papers:

  Janet, Coulomb, Chillet, Mas, "Simplified Magnetic Moment Method
  Applied to Current Transformer Modeling", IEEE Trans. Magn. 40(2),
  818-821, 2004.

  Janet, Coulomb, Chillet, Mas, "Magnetic Moment and Reluctance
  Network Mixed Method Applied to Transformer's Modeling", IEEE
  Trans. Magn. 41(5), 1428-1431, 2005.

## The target: current transformer

A current transformer (CT) has:
  - A primary winding carrying a possibly very large alternating
    current I_1 (kA range).
  - A small high-permeability iron toroidal core wrapped around the
    primary.
  - A secondary winding (many turns) burdened with a low-resistance
    "shunt" R_burden.  The induced secondary current I_2 = I_1 / N is
    measured across R_burden as a voltage signal.

The hard physics: high primary currents drive the core into deep
saturation; once the core saturates, large amounts of flux escape
into the air (LEAKAGE).  The output signal I_2(t) waveform is
distorted by both effects, and a designer needs to predict the
distortion as a function of geometry and core material.

## Magnetic-moment model applied to CTs

The core is meshed into n ferromagnetic elements with uniform
magnetization M_i in each (magnetization convention).  The total field
in element i is:

    H_i = H_0,i + sum_j f_ij M_j                                     (1)

where H_0,i is the external (Biot-Savart) field from the primary
winding evaluated at the center of element i, and f_ij is the
demagnetization tensor relating M_j to H_i.

Coupled with the BH law M_i = chi(H_i) H_i, solving (1) gives M_i.
Flux through the secondary coil is then a postprocessing integral
over the core volume.

Drawback: in CT geometry n must be large (100s to 1000s) to resolve
the long thin core; matrix assembly + solve becomes expensive.

## Simplified magnetic-moment model

Idea: use VERY FEW elements (3 to 5) for the entire core, with
magnetization direction PRESET in each element using symmetry.  In
the Janet "3-element model" the core is divided into 3 geometric
sectors, with the magnetization vector forced to lie along the median
plane of each sector (by symmetry of the toroidal CT).

The system becomes very small (3 unknowns), and computation reduces
to evaluating a 3x3 system at each saturation state.

Result: works WELL in deep saturation (error 3% vs FEM) -- because
when the core is heavily saturated, the directionality assumption is
correct and the magnetization is nearly uniform per sector.  But
FAILS at low induction -- because then the magnetization is small,
not uniform, and the coupling coefficients f_ij computed under the
saturation hypothesis are wrong.

## Simplified magnetic moment + Simplified RNA -- the Janet recipe

The fix: use RNA, which is accurate at low induction, to CORRECT the
magnetic-moment coupling coefficients f_ij.

1. Build a simplified RNA model of the same core (Fig. 6 Janet 2005):
   a few reluctances for the iron sectors + a few MMFs for the
   primary and secondary windings.

2. For a known LINEAR-regime operating point (small I_1), solve the
   RNA to get the mean magnetic field <H_i>_RNA and mean magnetization
   <M_i>_RNA in each element.

3. Compute the diagonal coefficients f_ii of the magnetic-moment coupling matrix
   by inverting the relation:

       <H_i>_RNA = H_0,i + f_ii <M_i>_RNA + sum_{j != i} f_ij <M_j>_RNA
                                                                     (2)

   where f_ij (j != i) are kept at their analytical FEM-computed
   values; only the SELF terms f_ii are corrected.

4. Optionally also correct the inductor matrix (the coefficients
   relating I_1 to H_0,i) by a scalar factor k chosen so that the
   magnetization in element 1 vanishes for a known zero-total-MMF
   condition (Janet 2005 eq. 7).

5. Use the corrected magnetic-moment system for the full saturation sweep.  Now
   both the low-induction and high-saturation regimes are captured by
   the same 3-5-element model.

## Performance (Janet 2005)

For their test CT:
  - Static error vs FEM: < 7% across the full operating range
    (compared to 15% for vanilla simplified magnetic-moment model and > 20% for pure 3-element
    magnetic-moment model without correction).
  - Static calculation time per point: 0.2 s (vs 1 min for full magnetic-moment model,
    20 min for FEM).
  - Response surface generation (140 x 140 grid): 35 min (vs 2 days
    for full magnetic-moment model, > 10 days for FEM).

The recipe trades modeling EFFORT (build both a small magnetic-moment model and a small
RNA, then perform the calibration) for run-time SPEED.  Worth it when
the model will be queried millions of times -- e.g. Monte Carlo
analysis, optimization, or as part of a real-time simulator.

## When to use RNA + magnetic-moment model vs full magnetic-moment model (Radia)

| Scenario                                  | Choice               |
|-------------------------------------------|----------------------|
| Need exact M(x) field map                 | HDiv-VIM / FEM     |
| Need fast secondary current waveform      | RNA + magnetic-moment model mixed      |
| Optimization over geometry parameters     | RNA + magnetic-moment model mixed      |
| Coupling to power-electronics SPICE model | RNA + magnetic-moment model mixed      |
| One-shot validation against measurement   | HDiv-VIM / FEM     |
| Sweep across 100+ load conditions         | RNA + magnetic-moment model mixed      |

## PEEC + magnetic-material coupling (Le-Duc 2013, COMPEL)

A different hybrid by Le-Duc, Chadebec et al. (G2Elab Grenoble):
combine PEEC for the conductors (driving AC currents) with a magnetic-material model for
the ferromagnetic material.  Both methods avoid meshing air, so the
hybrid retains the "mesh only what matters" virtue:

  - PEEC computes self / mutual partial inductances among conductor
    segments (Ruehli 1974).
  - The magnetic-material model provides the magnetic-material H field at each conductor
    segment, modifying the partial-inductance matrix.
  - The coupled system is solved in the frequency domain (time-
    harmonic) to give I_k and M_j simultaneously.

This is the natural counterpart of "RNA + magnetic-moment model": where RNA + magnetic-moment model
addresses CLOSED magnetic circuit + few-element 3D leakage, PEEC + magnetic-material model
addresses MULTI-CONDUCTOR ferromagnetic interaction (e.g. air-cored
inductor with a nearby iron yoke).

Cross-reference: `radia_mcp.peec.magnetic_policy` -- detailed PEEC + magnetic-material
formulation per Le-Duc 2013 with full derivation.
"""


ELECTROMECHANICAL_COUPLING = r"""
# Electromechanical coupling: RNA + ODE solver

When the magnetic system drives a mechanical motion (or vice versa),
the natural framework is to STATE-SPACE the coupled
electromagnetic-mechanical system and solve it with an ODE
integrator.

## Generic state-equation form

State vector: x = (z, v, I_1, I_2, ..., I_N) for N circuit branches.

Equations of motion + per-branch Kirchhoff voltage law (Lee 2005 eq.
1-2):

    m dv/dt = F_em(z, I) - m g - F_friction(v)                       (1)
    R_k I_k + d(lambda_k) / dt = V_k_applied                         (2)

with flux linkage lambda_k:

    lambda_k = sum_j L_kj(z) I_j                                     (3)

and electromagnetic force (Maxwell stress tensor OR virtual-work
gradient):

    F_em = 0.5 * sum_{i,j} I_i I_j d L_ij / dz                       (4)

The chain rule for d(lambda_k)/dt gives both an inductance-change and
a current-change contribution:

    d(lambda_k)/dt = sum_j L_kj(z) (dI_j/dt) + v * sum_j (dL_kj/dz) I_j
                                                                     (5)

The second term is the MOTIONAL EMF -- the back-emf induced by motion
of the conductor through a varying L(z).

## TEAM-28 levitation: the canonical demo

TEAM Workshop Problem 28 (Davey 2003, Lee 2005) is the standard
benchmark for the RNA + ODE pipeline:

| Quantity              | Value                                       |
|-----------------------|---------------------------------------------|
| Inside coil           | 960 turns, mean radius 39 mm, height 16 mm  |
| Outside coil          | 576 turns, mean radius 56 mm, height 16 mm  |
| Source                | I_src = 20 A peak, 50 Hz, series opposite   |
| Aluminum plate        | radius 53 mm, thickness 6.35 mm, mass 0.107 kg |
| Conductivity          | sigma = 3.4e7 S/m (Al)                      |
| Initial height        | 3.8 mm above coils                          |
| Steady levitated h.   | 11.3 mm (experimental)                      |

Expected response: at t = 0 turn on the AC, the plate jumps up,
overshoots, oscillates, and settles at ~11.3 mm.

The Lee 2005 procedure (covered in `lumped_extraction_fea` topic):
  1. Segment the plate into N_seg concentric rings.
  2. Sweep height z and FEM-extract L_kk(z), M_kj(z), and dL/dz.
  3. Time-march the state equations (RK4 explicit) with the per-step
     L(z), M(z), dL/dz looked up from cubic-spline interpolants.

Result (Lee 2005 Fig. 4-5, Table II):
  - 1 segment: 14.1% height error vs experiment
  - 2 segments: 5.5% error (matches the time-stepping FEM closely)
  - 5 segments: best agreement
  - Time: < 1 h vs 85 h for time-stepping FEM

The factor-of-100 speedup is the key result.  And note: changing
the input voltage / frequency / initial position does NOT require
a new FEM extraction -- the L(z), M(z) tables are reusable.

## Why RK4 (explicit) is enough here

The eddy-current time constants of the Al plate (a few ms) are
comparable to the mechanical time constants (10s of ms).  The system
is not stiff in practice.  A 4th-order Runge-Kutta with dt ~ 1e-5 s
is typical.

For STIFF systems (high-Q resonances, big disparity between EM and
mechanical timescales) one should switch to an implicit BDF or
Rosenbrock method.  Then the Jacobian of the L(z) terms must be
provided analytically or by finite differences.

## Coupling to a real circuit (SPICE-friendly)

The RNA + ODE solver framework dovetails directly with SPICE:
  - Treat the lumped magnetic device as a black box with port
    variables (V_port, I_port) and an internal state x.
  - Implement (V_port, I_port) -> (dx/dt, V_port_emf) as a SPICE
    behavioural source (Verilog-A or B-source).
  - The SPICE engine drives the device with the external power
    electronics; the device's internal state advances at the SPICE
    timestep.

This is how lab-scale Hane MEC + Cauer + Play models are typically
embedded in a switching-converter simulation -- the magnetic device
appears as a 1-port (transformer) or 2-port (coupled inductor) with
realistic iron loss + saturation + frequency-dependent eddy current.

## Cross-references

- `radia_mcp.motor.calc_motor_transient` -- general motor transient
  framework.
- `team28_reduced_model` topic -- TEAM-28 in further depth (modeling
  variations, comparison with BEM-FVM).
- public-safe curated corpus -- 10+ papers on TEAM-28
  variants (FVM, BEM, hybrid, POD, MOR).
"""


TEAM28_REDUCED_MODEL = r"""
# TEAM Workshop Problem 28: a deeper look

TEAM Workshop Problem 28 (Davey 2003) is the standard benchmark for
**eddy-current-driven electromechanical motion**: a thin aluminum
plate levitated above two coaxial coils.  It is small enough to be
solved by any method, complex enough to be a hard test of motion +
eddy-current + induction-heating coupling, and well-instrumented
experimentally.  Almost every electromechanical reduction method has
been tested against it.

## Problem statement (recap)

  - Two coaxial coils, inside / outside, wound with N_in = 960 and
    N_out = 576 turns, connected SERIES OPPOSITE (their fields nearly
    cancel along the axis -- a classical zero-net-flux design).
  - A 6.35 mm thick aluminum disk of mass 107 g sits 3.8 mm above
    the coils at t = 0.
  - AC current of 20 A peak at 50 Hz is applied.
  - The disk levitates, oscillates, and settles at h_steady = 11.3 mm
    (experimental observation, Davey 2003).

The challenge is reproducing the transient h(t) and v(t) curves with
acceptable accuracy.

## Reference solutions in public-safe curated corpus

| Method                              | Reference                  |
|-------------------------------------|----------------------------|
| Time-stepping FEM (full transient)  | Lee 2005 (reference value) |
| RNA + segmentation + RK4            | Lee 2005 (this topic)      |
| BEM-FEM coupled                     | "A modified BEM-FEM..."    |
| BEM-FEM with FMM acceleration       | "Fast multipole accel..."  |
| Meshless collocation                | "Application of meshless.."|
| 3D unstructured FVM                 | "Steady-State Analysis..." |
| Transient axisymmetric FVM          | "Transient Axisymmetric..."|
| POD-based MOR of eddy current       | "POD-based reduced-order.."|
| Matrix Interpolation MOR            | "Matrix Interpolation..."  |
| RL ladder with translational motion | "RL-Ladder Circuit..."     |

The Lee 2005 RNA + segmentation approach is the cheapest of these:
5 segments + RK4 ODE solver, 1 hour CPU vs 85 hours for the full
time-stepping FEM, 5% error vs experiment.

## What RNA + segmentation gets right and wrong

GETS RIGHT:
  - Steady-state levitation height (within 5%).
  - Damped oscillation period.
  - Per-segment eddy-current distribution if N_seg >= 5.
  - Sensitivity to input voltage / frequency / initial conditions
    (since the L, M tables are reusable).

GETS WRONG:
  - Initial transient overshoot can be inaccurate at N_seg = 1 (since
    the radial current distribution is poorly resolved).
  - Air friction not included by default -- final settling time is
    UNDERPREDICTED by ~10% (the experimental damping is partly
    aerodynamic, not electromagnetic).
  - Plate-tilt mode (rocking about a horizontal axis) is invisible to
    an axisymmetric formulation.  If the plate is not perfectly
    centered initially, it can tilt -- not modeled.

## Hybrid approaches in the literature

- **POD-based MOR** (Sato/Igarashi): Take snapshots from a few full
  FEM time integrations across the design space, extract POD modes,
  project the full FEM onto the POD subspace.  Result: 100x speedup
  vs full FEM but parametric design-space coverage requires N_modes
  snapshots.
- **Matrix Interpolation MOR** (Antoulas-style): Build a frequency-
  domain ROM at several geometry positions (z = z1, z2, ...) and
  interpolate the ROM matrices.  Fast at runtime, but ROM construction
  per z point is expensive.
- **RL ladder for translational motion**: Generalises Lee 2005 by
  expressing the position dependence of L(z), M(z) AS A CAUER LADDER
  in z (instead of cubic spline).  Useful for systems with many
  position-dependent eddy current loops.

## Take-home for RNA / MEC practitioners

1. If a small number of eddy-current loops dominate the motion (a few
   well-defined segments), Lee 2005 RNA + segmentation is hands-down
   the most efficient and easiest to maintain method.
2. If the geometry is so complex that you cannot identify the
   dominant eddy-current loops a-priori, switch to POD or Cauer
   ladder MOR (Kameari 2018) -- the modes are extracted automatically.
3. If you need TILT / PITCH / YAW dynamics, you cannot use an
   axisymmetric formulation and must go to full 3D FEM + 6-DOF rigid
   body dynamics.  TEAM-28 itself does not require this.

## Cross-references

- `electromechanical_coupling` topic -- general state-space recipe.
- `cauer_ladder_rna` topic -- CLN MOR for eddy-current fields.
- `radia_mcp.mor.snapshot.pod` -- POD MOR for parametric problems.
"""


TOPOLOGY_OPTIMIZATION = r"""
# RNA + Adjoint Variable Method for topology optimization

Topology Optimization (TO, since Bendsoe-Kikuchi 1988) lets the
optimizer pick not just dimensions but the actual SHAPE of the
magnetic material in a design region.  The classical approach uses
FEM as the inner solve + adjoint sensitivity, but each optimization
iteration costs a full FEM solve, so designs with 10^4 - 10^6 design
variables are slow.

Yin et al. (2023) replace the inner FEM with a fine-grid RNA, getting
a per-iteration cost reduction of ~10x.  The adjoint variable method
(AVM) is used to evaluate the gradient of the objective with respect
to ALL design variables at the cost of ONE EXTRA RNA solve.

Reference:

  Yin, Naidjate, Bracikowski, Pierquin, Trichet, "Topology
  Optimization of Magnetic Actuator based on Reluctance Network
  Modeling and Adjoint Variable Method", Compumag/CEFC, 2023.

## SIMP density formulation

Design variable: per-cell normalized density rho_i in [0, 1], where
rho_i = 1 means "iron" and rho_i = 0 means "air" in cell i.
Intermediate values mean a fictitious material with permeability

    mu(rho_i) = mu_0 (1 + (mu_r - 1) rho_i^n)                        (1)

(Yin eq. 5.)  The penalization exponent n > 1 discourages intermediate
densities by making them less effective per unit material -- standard
SIMP (Solid Isotropic Material with Penalization) trick.

The reluctance of cell i then depends on rho_i (via mu) and the optimizer
tunes rho to maximize the objective.

## RNA fine-grid model

The design region is meshed with a regular grid of (Nx x Ny) cells.
Each cell contributes 4 reluctances (Rl, Rr, Ru, Rd) and 2 flux
sources from any winding crossing the cell (eqs. 1-3 Yin).  The
global system:

    [Lambda(rho)] {eps} = {phi_src}                                  (2)

where Lambda is the reluctance matrix, eps the nodal magnetic
potential vector, phi_src the source vector.

## Objective: maximize armature force

The magnetic force on the armature is computed by Maxwell stress
tensor over its surface.  The optimization objective is:

    Maximize F(rho) = MST integral over armature surface
    subject to volume constraint sum rho_i = V_target

## Adjoint variable method (AVM) for sensitivity

Differentiate the state equation w.r.t. rho_i:

    d[Lambda]/drho_i * eps + [Lambda] * deps/drho_i = 0              (3)

The objective F can be considered an implicit function of (eps, rho):

    dF/drho_i = (dF/drho)_i + (dF/deps)^T deps/drho_i                (4)

Plug deps/drho_i = -[Lambda]^-1 (d[Lambda]/drho_i) eps into (4):

    dF/drho_i = (dF/drho)_i - (dF/deps)^T [Lambda]^-1 (d[Lambda]/drho_i) eps
              = (dF/drho)_i - lambda^T (d[Lambda]/drho_i) eps        (5)

with the ADJOINT VECTOR lambda defined by:

    [Lambda]^T lambda = dF/deps                                      (6)

Compute lambda ONCE (a single linear solve of the same size as the
forward RNA), then dF/drho_i for all i comes from (5) with the cost
of one matrix-vector product per i.

Cost comparison: finite-difference sensitivity costs N_design
extra RNA solves; AVM costs ONE extra solve.  For TO with N_design >
10^4, AVM is the only practical choice.

## Line search + filtering

Update rho in the steepest-descent direction with a line search to
ensure F monotonically increases.  After the update, apply a
density filter to ensure intermediate densities go to 0 or 1 cleanly
(feasibility-factor filter in Yin, alternatives include Heaviside
filter and smoothed projection).

## Result on the magnetic actuator example

Test geometry: C-shaped iron yoke + coil + movable armature.  Design
region: the air gap between the yoke poles.

  - Force evolution converges by ~iteration 8 (relative force
    improvement) and to a sharp 0/1 boundary by iteration ~35.
  - Optimal topology reduces leakage flux and concentrates flux into
    the armature -- the expected physical solution.
  - AVM speedup: orders of magnitude faster than finite-difference
    sensitivity for the same converged topology.

## When RNA-based TO is preferred over FEM-based TO

PREFER RNA-based TO:
  - Design region is naturally a 2D rectangular grid (PCB-style or
    laminated-sheet design).
  - Large number of design variables (10^4 - 10^6).
  - Objective is a lumped quantity (force, torque, flux linkage at a
    port) rather than a local field distribution.
  - Linear or weakly-nonlinear material.

PREFER FEM-based TO:
  - Curved boundaries or complex 3D shapes.
  - Strongly-nonlinear / saturating iron.
  - Objective involves local stress or local heating.
  - High accuracy required (5% level rather than topology-only).

## Cross-references

- `mec_basics` topic -- grid-based reluctance assembly.
- `reluctance_network_construction` topic -- Yin grid-cell formula.
"""


DYNAMIC_HYSTERESIS = r"""
# Dynamic hysteresis in MEC -- Play + Cauer (Sugahara-Nakamura-Hane lineage)

Static-BH MEC (Hopkinson with mu(B)) is fine for slowly-varying
fields, but FAILS at PWM-driven hysteresis loops with minor-loop
patches at carrier frequencies of 1-10 kHz.  The eddy-current
distribution inside the lamination is not uniform at those
frequencies (skin depth < lamination thickness), and the apparent
permeability becomes complex and frequency dependent.

The Sugahara-Nakamura-Hane lineage solves this with a dynamic
hysteresis MEC that combines:
  - PLAY MODEL for the DC hysteresis (Bobbio-Visone 1997, extended to
    Tohoku-Kindai shape-function decomposition).
  - CAUER LADDER NETWORK (Cauer-I) for the eddy-current skin effect
    inside the lamination.
  - Separate dependent source for ANOMALOUS eddy-current loss
    (Bertotti excess loss).

Key references:

  Hane, Nakamura, "Dynamic Hysteresis Modeling for Magnetic Circuit
  Analysis by Incorporating Play Model and Cauer's Equivalent Circuit
  Theory", IEEE Trans. Magn., DOI 10.1109/TMAG.2020.3004355, 2020.

  Tanaka, Nakamura, Ichinokura, "Magnetic Circuit Model combined with
  Play Model Obtained from Landau-Lifshitz-Gilbert Equation", JPCS
  903, 012047, 2017.

  Shindo, Noro, "Simple Circuit Simulation Models for the Eddy Current
  in Magnetic Sheets and Wires", IEEJ Trans. FM 134(4), 173-181,
  2014.

  Shindo, Miyazaki, Matsuo, "Cauer Circuit Representation of the
  Homogenized Eddy-Current Field Based on the Legendre Expansion for
  a Magnetic Sheet", IEEE Trans. Magn. 52(3), 6300504, 2016.

## The MEC per-branch structure

Each branch in the RNA is replaced by:

    Branch i = [ Play hysterons ] in series with [ Cauer ladder ]

The PLAY HYSTERON block (Hane 2020 Fig. 1) takes the branch flux Phi
as input and produces the MMF F as output via:

    H = sum_k zeta_k(eta_k, B)                                       (1)

where eta_k are play widths (Tesla), zeta_k are play hysterons, and
the sum is weighted by shape functions f_k(r) derived from a family
of MEASURED DC hysteresis loops at different B_max amplitudes
(Hane Fig. 2).

The CAUER LADDER block (Hane 2020 Fig. 6, 7) is a Cauer-I network
truncated to 2 or 3 stages with parameters (Hane eq. 2-3):

    L_2n+1 = mu * geometric_factor                                   (2)
    R_2n   = 4 / (sigma d^2)                                         (3)

where d is lamination thickness, sigma conductivity, and mu the
LINEARIZED differential permeability evaluated on the centerline of
the DC hysteresis loop at B_m = 1 T.  This linearization is a
practical recipe Hane validates empirically on grain- and non-
oriented Si-Fe.

The CLASSICAL eddy-current loss equals 1/L_1 in the classical limit
(first stage of ladder = single inductance).  Adding more stages
captures the skin-effect-driven reduction at high frequency.

The ANOMALOUS eddy-current loss is represented as a separate
dependent source proportional to (dB/dt)^1.5 with coefficient gamma_2
calibrated from low-frequency core-loss measurements (Hane Fig. 3 and
Table I).

## Time stepping

At each time step t:

  1. Read winding currents from external SPICE circuit -> MMF sources
     on the source branches.
  2. Solve nonlinear KCL at all nodes (modified-Newton with current
     branch tangents).
  3. For each branch:
     - Update the Play state (which hysterons are pushed)
     - Update the Cauer state variables (currents in each ladder
       inductor)
     - Sum the resulting H, the dB/dt contribution, and the anomalous
       term.
  4. Compute total flux per branch + per-branch iron-loss breakdown:
       P_total = P_hyst + P_classical + P_anomalous

The Play state update is O(K) per branch (K = number of hysterons,
typically 5-20).  The Cauer state update is O(N_stage) per branch
(N_stage = 2-5).  Both are constant-cost-per-step, no matrix solve
needed for the hysteresis-loss part itself.

## Iron-loss decomposition (Hane Fig. 9)

The output of the model gives:

  - P_hyst (Joule loss in static hysteresis loop) -- frequency
    independent.
  - P_classical (Joule loss in distributed eddy currents) -- scales
    as f^2 below the skin-depth corner, falls off above.
  - P_anomalous (Joule loss in domain-wall fluctuations) -- scales
    as f^1.5.

This decomposition is the SUGAHARA-NAKAMURA-HANE selling point: not
just total loss but a physically-interpretable breakdown that
guides material selection (high-Si steel reduces P_classical;
domain-refined steel reduces P_anomalous).

## Complex permeability output (Hane Fig. 10)

By exciting the model with a small sinusoidal H at frequency f and
fitting the response, one extracts:

    mu*(f) = mu'(f) - j mu"(f)                                       (4)

The model reproduces measured mu'(f) and mu"(f) across 10 Hz - 10 kHz
on grain- and non-oriented Si-Fe.  This validates the model not just
for total loss but for the full frequency-domain magnetic response,
which is what a power-electronics circuit simulator actually sees.

## Forward extension: full motor RNA

The Hane "1912 マグ研" lab paper applies the per-branch Play + Cauer
structure to a multi-branch RNA of a permanent-magnet motor / three-
phase variable inductor.  Each branch's hysteresis and eddy-current
state evolves independently; the nonlinear KCL couples them through
the flux-conservation constraint.

This is THE most physically detailed RNA available for PWM iron-loss
prediction in rotating machines, and it is the production model of
the Tohoku-Kindai collaboration.

## Cross-references

- `cauer_ladder_rna` topic -- the Cauer ladder by itself.
- `radia_mcp.magnetic_materials.hysteresis_models.play` -- Play model
  forward / inverse evaluation API.
- `radia_mcp.magnetic_materials.hysteresis_models.lab_core` -- the
  Tohoku-Kindai energy / play core lineage.
- public-safe curated corpus -- the
  seminal lab paper.
- public-safe curated corpus --
  English IEEE version.
"""


VS_PEC_PEEC = r"""
# Terminology: RNA vs PEC vs PEEC vs MEC vs MoM vs HDiv-VIM

These acronyms are often confused.  Here is the canonical split as
used in the Sugahara-Kindai lab and in the literature cited above:

## RNA -- Reluctance Network Analysis

A magnetic-circuit method.  Discretize a magnetic structure into
LUMPED reluctances + MMF sources + flux sources, solve via Kirchhoff
laws (KCL or KVL).  Domain: magnetostatics + slow-frequency
quasi-static.  Output: branch fluxes, port inductance, force on
armature, iron loss with optional Play / Cauer extensions.

Synonyms: MEC (Magnetic Equivalent Circuit), magnetic-circuit method,
permeance network.

## MEC -- Magnetic Equivalent Circuit

Identical to RNA.  Different communities use different names:
  - "RNA" preferred in Japanese / Tohoku-Kindai literature.
  - "MEC" preferred in American IEEE Energy Conv. literature (Slemon,
    Ostovic, Sudhoff, Pekarek, Derbas).
  - "Permeance network" preferred in classical Roters 1941 framing.

In radia_mcp the topic name `rna_mec` captures both.

## PEC -- Partial Element Equivalent Circuit (DEPRECATED ALIAS)

This acronym is sometimes used as shorthand for PEEC; rarely seen in
serious literature.  Prefer PEEC.

## PEEC -- Partial Element Equivalent Circuit

An ELECTRICAL-circuit method by Ruehli (1974) for solving Maxwell
equations in the time-harmonic / quasi-static regime.  Discretize
conductors into segments + panels; compute partial inductances L_ij,
partial capacitances P_ij, and segment resistances R_i analytically
from segment geometry.  The full Maxwell solve becomes an MNA
(modified nodal admittance) solve of an EQUIVALENT CIRCUIT made of
(R, L, P, M) elements.

Domain: conductor systems with AC currents, full 3D, frequency-
dependent skin effect via SIBC.

In radia_mcp: see `radia_mcp.peec`.

## MoM -- Method of Moments

A general numerical method (Harrington 1968) for solving integral
equations.  In electromagnetics, "MoM" usually means an EFIE / MFIE
solver for SCATTERING / radiation problems (full-wave, Helmholtz
kernel).

Sometimes "MoM" is used loosely to mean any integral-equation method,
including PEEC.  Strictly, MoM = Galerkin / collocation
discretization of an integral operator.

## HDiv-VIM -- Radia magnetic-material route

Radia's current magnetic-material route is HDiv-VIM: express the
magnetic material problem in NGSolve-compatible HDiv language, build
the charge-Gram interaction, and keep the route compatible with
reduced FEM and higher-order curved meshes.

In radia_mcp: see `radia_mcp.radia_ngsolve.hdiv_vim`.

## Cross-method comparison table

| Acronym | Unknowns         | Domain               | Material        |
|---------|------------------|----------------------|-----------------|
| RNA/MEC | branch flux phi  | magnetostatics + QS  | iron (BH curve) |
| PEEC    | seg current I    | air + conductors, QS | linear (or SIBC)|
| MoM     | surface curr J   | full-wave            | PEC, lossy diel.|
| HDiv-VIM| HDiv / charge DOF| magnetostatics       | nonlin iron     |
| FEM     | nodal A or H     | any domain           | any             |
| BEM     | surface phi or J | any domain           | any             |

## Which to use when

- Pure soft iron + PM, no eddy currents -> **HDiv-VIM** (Radia).
- Conductor inductance extraction, possibly with skin effect ->
  **PEEC** (Radia).
- 1-D laminated sheet eddy current -> **Cauer ladder** (Hane 2020).
- Quick electric-machine concept design with saturation -> **RNA / MEC**.
- Force on a moving conductor with eddy currents -> **RNA + segmentation +
  RK4** (Lee 2005 TEAM-28).
- Topology optimization of an actuator -> **grid RNA + AVM** (Yin 2023).
- PWM iron-loss prediction in laminated machines -> **dynamic
  hysteresis MEC** = RNA + Play + Cauer (Hane / Sugahara-Nakamura
  lineage).
- High frequency (1 MHz+) scattering or antennas -> **MoM** (not
  Radia, use a dedicated solver).
- Detailed local field map in complex 3D -> **FEM** (NGSolve).

The Sugahara-Kindai-Tohoku research program covers RNA + HDiv-VIM +
PEEC + Cauer + Play in a unified framework; for any specific problem
the right combination is usually obvious from the table above.
"""


def get_knowledge(topic: str = "overview") -> str:
    """Dispatch RNA / MEC topics.

    Topics (see TOPICS dict for one-line summaries):
        overview                       (DEFAULT)
        mec_basics
        nodal_vs_mesh_analysis
        reluctance_network_construction
        lumped_extraction_fea
        cauer_ladder_rna
        rna_magnetic_coupling
        electromechanical_coupling
        team28_reduced_model
        topology_optimization
        dynamic_hysteresis
        vs_pec_peec
        all
    """
    topic = topic.lower().strip()
    aliases = {
        "intro": "overview",
        "mec": "mec_basics",
        "basics": "mec_basics",
        "reluctance": "mec_basics",
        "nodal_vs_mesh": "nodal_vs_mesh_analysis",
        "nodal_mesh": "nodal_vs_mesh_analysis",
        "derbas": "nodal_vs_mesh_analysis",
        "construction": "reluctance_network_construction",
        "network_construction": "reluctance_network_construction",
        "fea_extraction": "lumped_extraction_fea",
        "lumped_extraction": "lumped_extraction_fea",
        "lee": "lumped_extraction_fea",
        "cauer": "cauer_ladder_rna",
        "cln": "cauer_ladder_rna",
        "kameari": "cauer_ladder_rna",
        "rna_magnetic": "rna_magnetic_coupling",
        "transformer": "rna_magnetic_coupling",
        "janet": "rna_magnetic_coupling",
        "electromechanical": "electromechanical_coupling",
        "team28": "team28_reduced_model",
        "team_28": "team28_reduced_model",
        "topopt": "topology_optimization",
        "topology_opt": "topology_optimization",
        "to": "topology_optimization",
        "yin": "topology_optimization",
        "dynamic": "dynamic_hysteresis",
        "play_cauer": "dynamic_hysteresis",
        "hysteresis_mec": "dynamic_hysteresis",
        "hane": "dynamic_hysteresis",
        "terminology": "vs_pec_peec",
        "acronyms": "vs_pec_peec",
        "vs_peec": "vs_pec_peec",
    }
    topic = aliases.get(topic, topic)

    sections = {
        "overview": OVERVIEW,
        "mec_basics": MEC_BASICS,
        "nodal_vs_mesh_analysis": NODAL_VS_MESH_ANALYSIS,
        "reluctance_network_construction": RELUCTANCE_NETWORK_CONSTRUCTION,
        "lumped_extraction_fea": LUMPED_EXTRACTION_FEA,
        "cauer_ladder_rna": CAUER_LADDER_RNA,
        "rna_magnetic_coupling": RNA_MAGNETIC_COUPLING,
        "electromechanical_coupling": ELECTROMECHANICAL_COUPLING,
        "team28_reduced_model": TEAM28_REDUCED_MODEL,
        "topology_optimization": TOPOLOGY_OPTIMIZATION,
        "dynamic_hysteresis": DYNAMIC_HYSTERESIS,
        "vs_pec_peec": VS_PEC_PEEC,
    }
    if topic in sections:
        return sections[topic]
    if topic == "all":
        return "\n\n".join(sections[k] for k in [
            "overview",
            "mec_basics",
            "nodal_vs_mesh_analysis",
            "reluctance_network_construction",
            "lumped_extraction_fea",
            "cauer_ladder_rna",
            "rna_magnetic_coupling",
            "electromechanical_coupling",
            "team28_reduced_model",
            "topology_optimization",
            "dynamic_hysteresis",
            "vs_pec_peec",
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        + ", ".join(sorted(sections.keys()))
        + ", all."
    )
