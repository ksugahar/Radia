# High-Order HCurl Response Compression for CLN/VIM Eddy-Current Models

This note records the theory direction for the hybrid eddy-current model:
use a high-order NGSolve `HCurl` space as the parent space, keep only the
external-field-visible response subspace, and assemble a passive reduced
VIM/CLN model enriched by surface-Omega/SIBC branches.

The intended presentation angle is theoretical.  Numerical experiments should
illustrate the construction, not carry the main claim.

---

## 1. Core Claim

Let `V_p subset H(curl, Omega_c)` be a high-order conductor space and let
`B` be the small set of low-order external-field ports.  The useful eddy-current
degrees of freedom are not all DoFs of `V_p`, but the response-visible Krylov
subspace

```text
R_{p,n}(B) =
span{ K_p^{-1} B,
      K_p^{-1} M_p K_p^{-1} B,
      ...,
      (K_p^{-1} M_p)^{n-1} K_p^{-1} B }.
```

Here `K_p` is the magnetic stiffness / inductive operator and `M_p` is the
conductive mass / resistance-side operator in the high-order `HCurl` space.
The integer `p` controls how accurately the parent space can represent spatial
skin-effect modes.  The integer `n` controls how many CLN/Cauer moments are
kept.  They are related, but they are not the same parameter.

The reduction is therefore:

```text
high-order HCurl parent space
    -> response-visible subspace
    -> topology guard for eddy-current loop closure
    -> curl(T) sampled current basis
    -> passive reduced VIM / CLN circuit
```

Terminology used in this project:

| Term | Meaning |
|---|---|
| **Eddy-Visible Response Space (EVRS)** | The retained response subspace visible through the external-field ports and the chosen CLN/Krylov depth. |
| **Eddy-Invisible DoF (EID)** | A high-order eddy-current DoF or direction that exists in the parent `HCurl(p)` space but is not visible to the selected low-order external-field ports within the retained response order. |
| **eddy bubble** | Informal short name for an EID.  This is not the same as an NGSolve element bubble / `LOCAL_DOF`. |

This is stronger than using a low-order edge-element model directly.  With
`p=1`, later Cauer stages can be limited by the parent space before the reduced
model has a chance to express the physics.  With `p=3,4,...`, the parent space
contains sharper eddy-current distributions, and the response compression can
discard the EIDs after they have done their job.

The compression must not be purely local.  Eddy currents flow in closed loops,
so a basis direction that looks small on one element can still be the bridge
that closes a loop through a neighboring conductive element.  The first
production rule is therefore face-topological:

| Face adjacency | Basis consequence |
|---|---|
| conductor-air / exterior conductor face | allow surface-Omega/SIBC boundary basis and no normal current leakage |
| conductor-conductor face with the same conductive material | preserve loop-bridge directions across the face |
| conductor-conductor face with different conductive materials | preserve loop-bridge directions and material-interface continuity |
| conductor-insulator face | treat as an exterior current-blocking face, but do not call it an SIBC half-space unless the physics supports that approximation |

Radia records this distinction with `ClassifyNgsolveEddyTopology`.  Air-touching
conductor faces, classified as `conductor-air` or `conductor-exterior`, are the
SIBC faces.  `NgsolveEddyDofPolicy` then maps those face roles to HCurl DoF
masks:

```text
conductor-air/exterior DoFs  -> SIBC/surface-Omega reduction candidate
conductor-insulator DoFs     -> non-SIBC boundary trace candidate
conductor-conductor          -> topology-protected loop bridge
remaining interior           -> ordinary EVRS / eddy-bubble reduction candidate
```

This classification is face-local; it is not permission to delete every
element that does not touch air.  At low frequency the skin depth may be
comparable with the conductor thickness, so the current fills the volume.
Conductor-conductor adjacency also carries the cycle topology that closes the
current.  The removable object is therefore an interior response direction
that is invisible to the selected ports and observables, not an interior
element merely because it is remote from the exterior boundary.

The production reduction is performed by class, not by a single global cutoff:

```text
air/exterior surface class:  replace by surface-Omega/SIBC basis or boundary ladder
non-air trace class:         keep as a boundary trace condition, not SIBC
bridge class:                reduce conductor-conductor face bridges by a graph cycle basis
interior class:              compress aggressively by EVRS / CLN Krylov depth
```

Validation JSON stores both `topology_diagnostics` and
`dof_policy_diagnostics`.  `EddyConductorGraph` records the conductive-element
adjacency graph and its fundamental cycle basis.  `EddyReductionPlan` records
the class-wise production plan so later EVRS/HACApK/BEM runs can show how many
faces were treated as SIBC candidates, how many conductor-conductor bridge DoFs
were replaced by cycle-basis modes, and how many DoFs remain ordinary EVRS
reduction candidates.

`EddyBubbleDecomposition` is the named production object for this split.  It
exposes the structurally retained mask

```text
structural_keep = SIBC surface + non-SIBC trace + bridge/cycle
```

and the removable bulk mask

```text
eddy_bubble_candidate = ordinary interior bulk EVRS candidate.
```

This is the point where "eddy bubble" becomes an implementation contract rather
than a loose nickname.

`NgsolveBridgeCycleCurrentBasis` then turns those graph cycles into a coarse
VIM-compatible current basis.  It samples each conductor-conductor face at its
face center, orients the current between adjacent element centroids, and uses
the dual-volume weight `area(face) * centroid_distance`.  This is a bridge
between topology and physics: the cycle rank no longer remains only a count,
but becomes an actual reduced current basis that can enter `AssembleHybridVIM`.

This split is now also recorded as a Mathematica certificate in
`packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/eddy_topology_reduction.wls`.
The symbolic statement is: the face role is determined by the two neighboring
material labels; only conductor-air/exterior faces are SIBC/surface-Omega
traces; conductor-insulator faces are non-SIBC boundary traces;
conductor-conductor faces form a conductive adjacency graph; if `B` is its
oriented incidence matrix and `Z` is a cycle basis, then `B Z = 0`, so
`j_bridge = Z gamma` preserves cell-wise current continuity exactly.  Thus
`p` should be argued as a parent-space admissibility order, not selected from a
desktop p-sweep alone.

The companion order ledger
`packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/eddy_parent_order_admissibility.wls`
makes that statement explicit:

```text
p >= max(r_bulk, r_bridge, r_surface),
n_reduced = n_EVRS + cycle_rank * dim P_r_bridge(face) + n_surface.
```

Here `r_bulk` is the spatial current degree that the bulk EVRS/CLN moment
construction needs, `r_bridge` is the polynomial trace degree retained on
conductor-conductor bridge faces, and `r_surface` is the trace degree retained
on SIBC/surface-Omega faces.  A surface-Omega scalar potential of degree
`r_surface+1` generates a surface current of degree `r_surface`.  Therefore
`p=6` is a symbolic claim only when one of these requirements reaches 6; in the
current coarse bridge-cycle smoke it is better described as admissible and
conservative, not optimal.

The same ledger is available from Python as `EddyParentOrderLedger`.  This is
the production-side guard against accidentally hard-coding `p=6` as a rule.  A
typical diagnostic is:

```python
ledger = vim.EddyParentOrderLedger(
    bulk_degree=4,
    bridge_trace_degree=0,
    surface_current_degree=2,
)
ledger.diagnostics(parent_order=6, evrs_rank=24, cycle_rank=7, surface_modes=3)
```

For that ledger the required parent order is 4, so `p=6` is admissible with an
excess of 2.  If a higher-order bridge trace is needed, the mode count changes
as `cycle_rank * EddyTracePolynomialDim(r_bridge)` and the required parent
order rises only when `r_bridge` is the dominant requirement.

---

## 2. Coexistence with HDiv-VIM

The construction should coexist with Radia's existing HDiv-VIM path, in the same
spirit as an ELF-style mixed formulation.

The reason is the de Rham map:

```text
T in H(curl)  --curl-->  J = curl T in H(div),     div J = 0.
```

The eddy-current branch is generated in `HCurl`, but the physical current basis
used by VIM is the sampled `curl(T)`, which is an HDiv-compatible solenoidal
current.  This makes the method a neighbor of the existing HDiv-VIM machinery,
not a competing vocabulary.

The symbolic Mathematica check
`packages/radia-mcp/src/radia_mcp/mathematica/basis_functions/evrs_derham_bridge.wls`
locks the local FEEC identity used here: an arbitrary retained EVRS combination
inside the HCurl parent space still satisfies `div(curl T_EVRS)=0`, and the
surface-Omega branch `K = n x grad_Gamma Omega` is tangential and
surface-solenoidal on a reference flat face.

In that Mathematica derivation the eddy variable is fixed as `T`.  On the
Whitney tetrahedron it constructs the discrete chain

```text
phi_coeff  --G-->  T_coeff  --C-->  J_coeff  --D-->  rho_coeff
```

and verifies

```text
C G = 0,          D C = 0,
T_coeff = Q z,   J_coeff = C Q z,   D J_coeff = 0.
```

It also verifies the T-method matrix form

```text
R_T = C^T M_R C,       L_T = C^T M_L C
```

for generic symmetric current-space matrices `M_R`, `M_L`: the matrices are
symmetric, annihilate gradient gauge directions, and after EVRS projection obey

```text
(Q z)^T R_T (Q z) = (C Q z)^T M_R (C Q z).
```

The current-space port RHS has the same gauge invariance,
`G^T C^T P = 0`, and the reduced port RHS is `Q^T C^T P`.

The practical split is:

| Role | Space / model |
|---|---|
| Magnetization, flux, soft iron, open-region demag | existing Radia HDiv-VIM |
| High-order eddy-current response generation | NGSolve `HCurl` parent space |
| Reduced volume current basis | `J_i = curl T_i`, HDiv-compatible |
| Surface skin branch | surface-Omega basis `K_i = n x grad_Gamma Omega_i` |
| Exterior interaction | VIM Laplace single-layer backend; BEM/H-matrix can be used only as projection accelerators |

This gives a natural rotor-stator path: Radia/HDiv-VIM can own magnetic bodies
and open-region kernels, while the stator or conductor response is compressed
from high-order NGSolve `HCurl` spaces.

---

## 3. Relation to CLN Stages

For a linear conductor, the port transfer function obtained from a finite
element discretization is a positive-real rational function.  A Cauer ladder is
a continued-fraction realization of the same positive-real transfer function.

The response basis above is a block Krylov basis.  In exact arithmetic, a
Galerkin projection onto the first `n` Krylov blocks matches the corresponding
low-frequency moments of the transfer function.  A Lanczos/Stieltjes
continued-fraction realization of those moments gives the Cauer ladder.  This
is the algebraic bridge between

```text
high-order FE matrices  <->  response Krylov basis  <->  CLN stages.
```

The parent response basis must then be quotiented by the current map.  Let
`V` contain the HCurl Krylov vectors, `C` sample `curl(T)`, and `W` contain the
volume quadrature weights.  Radia forms

```text
G_J = (C V)^* W (C V).
```

Directions in the relative null space of `G_J` are independent parent-T
coordinates but carry no independent eddy current, loss, or VIM interaction.
`CompressHCurlResponseInCurrentGram` removes them and applies the whitening
transform to both `V` and `C V`.  This second quotient is the current-space
Eddy Bubble removal; topology-aware bridge/SIBC protection is applied after it.

The key design rule is therefore not "use p=4" or "use four stages".  The rule
is:

```text
Choose p high enough that the first n CLN stages are stable under p-refinement.
Then choose n by the frequency band and by the desired circuit accuracy.
```

A useful research diagnostic is the table

```text
p = 1,2,3,4,5,...
n = 1,2,3,4,5,...
check: positivity of R_k,L_k; p-refinement of moments; port transfer error.
```

The first runnable Radia validation lane for this table is
`validation_test/cln/evrs_pn_convergence.py`.  It builds a small unit-box
`HCurl(p)` parent problem, compresses to EVRS depth `n`, samples `curl(T)`,
assembles the reduced VIM, and writes `evrs_pn_convergence_smoke.json` with
DoF compression and port-admittance error relative to the highest `(p,n)` case
in the run.  This is a desktop smoke, not benchmark evidence; non-monotone rows
are useful diagnostics of sampling order, kernel regularization, or an
insufficient reference case.

The companion lane `validation_test/cln/evrs_sibc_mixed_schur.py` adds
surface-Omega/SIBC modes to the same bulk EVRS basis and checks the mixed
Galerkin Schur complement.  It accepts `--orders 4,5,6` so the same validation
can compare candidate parent orders against the highest-p reference case.  In
the current p=6 smoke, the 1226 active parent
DoFs split into 258 SIBC-surface DoFs, 386 loop-bridge DoFs, and 600 ordinary
EVRS candidates.  The conductor graph has cycle rank 7, so the class-wise
tri-block system at the reference depth is 24 bulk EVRS coordinates, 7
bridge-cycle coordinates, and 3 surface coordinates: 34 modes in total.  The
n=11 bulk basis is already within about 0.11% in the mixed port admittance,
while the Schur-reconstructed mixed-system residual is about 1e-13.

The question is not whether `p=6` is inherently needed.  The production choice
should be the smallest parent order that satisfies the symbolic order ledger
and leaves the retained observables invariant under further p-refinement.  For
many smooth low-frequency motor cases that may be `p=4` or `p=5`; `p=6` is a
conservative parent-space smoke point unless `r_bulk`, `r_bridge`, or
`r_surface` actually reaches 6.  A thinner skin depth, sharper corner field, or
higher-order SIBC correction can raise the ledger requirement.

The first `p=4/5/6` smoke supports that interpretation.  With the default
ledger (`required_parent_order=4`) and a smooth unit-box conductor,
`validation_test/cln/evrs_sibc_p456_depth20_smoke.json` gives `p=4,n=18`
within 0.25% of the `p=6,n=20` reference, while the active parent DoFs are 428
instead of 1226.  Thus `p=6` is not justified as a universal default by this
case.  It remains a useful conservative reference while sharper geometry,
SIBC-trace enrichment, and motor-grade meshes are being tested.

The first sharper-geometry smoke points in the opposite direction.  With
`--geometry notched-box`, a re-entrant exterior corner, and a one-frequency
comparison against `p=6,n=22`, the p=4 and p=5 rows stay about 2--3% away from
the p=6 reference even at the same EVRS depth.  This is not yet a theorem, but
it is the first numerical evidence that the conservative `p=6` parent space can
matter once the field has corner-driven spatial content.  The practical rule is
therefore: use the ledger for admissibility, then run p-refinement on the
retained port/SIBC observables for the actual geometry class.

The field-level comparison is stronger.  In
`validation_test/cln/evrs_current_field_compare.py`, the full p=6 HCurl parent
solve is used as the reference and the sampled current `J = curl(T)` is compared
directly.  On the same notched-box geometry, p=6 can be reduced from 3557 active
parent DoFs to 44 EVRS coordinates while keeping the current L2 error at about
`1.6e-6`, Joule-loss error below `1e-7`, and corner-local current error around
`1e-6`.  A p=4 parent reduced to the same 44 coordinates is still about 21% away
in current L2 and about 15% away in Joule loss.  Thus the high-order p=6 parent
space is doing real work, but only a tiny EVRS subspace of it has to survive in
the production eddy-current basis.

---

## 4. Passive VIM Projection

Let the reduced current basis be stored in the columns of `Q`.  The reduced
matrices are

```text
R_r = Q^* R Q,       L_r = Q^* L Q.
```

If `R` is Hermitian positive semidefinite and `L` is Hermitian positive
semidefinite, then the projected matrices have the same property.  Therefore

```text
Z_r(s) = R_r + s L_r
```

is passive for `Re(s) >= 0`.  With a surface branch,

```text
Z_r(s) = R_r + s L_r + Z_s(s) M_Gamma
```

is still passive if `Re(Z_s(s)) >= 0` and `M_Gamma` is Hermitian positive
semidefinite.

This is the reason the reduced model should be assembled as a Galerkin/VIM
projection, not as an ad hoc fitted circuit.  Positivity and reciprocity are
structural.

The production assembly order is important.  A dense p=6 parent-space Green
Gram is not formed.  NGSolve first assembles the sparse high-order HCurl
stiffness, mass, and port actions.  Static condensation and EVRS then produce
the retained `J = curl(T)` modes, bridge-cycle modes, and SIBC modes.  Only this
small retained current basis enters the VIM Gram.  Near interactions require
singularity-aware or adaptive integration; far interactions may use the same
Laplace-kernel ACA/H-matrix backend as the other Radia VIM operators.  Thus p=6
controls parent-space resolution without making the dense Gram dimension p=6.

---

## 5. SIBC Bridge

The SIBC branch should be viewed as a boundary Cauer branch.

For a local half-space conductor,

```text
Z_s(s) = sqrt(mu s / sigma)
```

is a positive-real surface impedance for `Re(s) >= 0`.  For a finite-thickness
slab, the surface impedance is a hyperbolic-function variant of the same
diffusion solution.  Rational Cauer/Pade approximants of these one-dimensional
diffusion impedances give passive ladder branches.

A surface-Omega basis represents the tangential surface current as

```text
K = n x grad_Gamma Omega.
```

The zeroth-order surface branch uses only the surface mass matrix
`M_Gamma`.  Higher-order SIBC corrections introduce geometry and surface
diffusion operators, such as curvature terms and a Laplace-Beltrami stiffness.
Thus the bridge is:

```text
local diffusion SIBC
    -> positive-real boundary impedance Z_s(s)
    -> Cauer/Pade boundary ladder
    -> surface-Omega reduced basis
    -> same VIM/BEM exterior interaction matrix.
```

This does not imply that "SIBC order = CLN stage".  SIBC order is an asymptotic
order in skin depth and curvature.  CLN stage is a rational approximation order
for a transfer function.  They meet when the SIBC operator is itself realized as
a passive boundary ladder.

For a linear conductor, an exact discrete DtN map can eliminate volume
unknowns and retain only boundary traces.  The resulting operator is dense and
frequency dependent; a local SIBC is its boundary-layer approximation, valid
when penetration is small relative to the local thickness and curvature scale.
Only conductor-air/exterior faces are candidates for that SIBC replacement.
Conductor-conductor faces remain topology-bearing interfaces.

Nonlinearity does not make a surface model mathematically impossible, but it
changes the contract.  A saturating or hysteretic conductor has a field- and
history-dependent DtN map, so one precomputed linear boundary operator is no
longer sufficient.  A nonlinear ESIM/SIBC cell model remains useful as an
explicit reduced approximation.  The robust reference and the production path
when penetration, corner fields, or material state matters is the volume VIM
branch, with response reduction performed after the material operator has been
assembled or updated.

---

## 6. Implementation Status in Radia

The initial API is in `src/radia/vim/_eddy_hybrid.py`.

Current primitives:

| API | Role |
|---|---|
| `EVRSBasis` | named retained Eddy-Visible Response Space with EID/eddy-bubble diagnostics |
| `CompressHCurlResponseInCurrentGram` | removes response directions in the relative null space of sampled `curl(T)` and current-Gram orthonormalizes the retained parent/current basis |
| `NgsolveBlockKrylovBasis` | NGSolve matrix/vector bridge into response-visible basis |
| `NgsolveStaticCondensedBlockKrylovBasis` | NGSolve `LOCAL_DOF` static condensation followed by EVRS/EID response compression |
| `NgsolveHCurlCurlBasis` | maps response coefficients to sampled `curl(T)` current modes |
| `NgsolveHDivMagnetizationBasis` | maps HDiv coefficient vectors to sampled magnetization modes |
| `NgsolveSurfaceOmegaBasis` | maps surface gradients to `K = n x grad_Gamma Omega` |
| `ClassifyNgsolveEddyTopology` | classifies conductor-air/exterior faces versus conductor-conductor loop bridges before pruning |
| `EddyConductorGraph` | builds the conductive-element graph and fundamental cycle basis for bridge reduction |
| `NgsolveEddyDofPolicy` | maps topology roles onto HCurl DoF masks: SIBC surface, non-SIBC trace, loop bridge, and ordinary bulk EVRS candidate |
| `EddyReductionPlan` | records the class-wise production reduction plan and optional reduced-mode estimate |
| `EddyBubbleDecomposition` | named production split for topology-aware eddy bubbling: structural keep versus bulk eddy-bubble candidates |
| `EddyBubbleReduction` / `NgsolveEddyBubbleReduction` | construct the eddy-bubbling split from an existing policy or directly from an NGSolve HCurl space |
| `EddyBubbleHCurlBasis` | production basis object: EVRS coefficients, sampled `J=curl(T)`, eddy-bubbling diagnostics, VIM assembly, and HDiv-MMM coupling hook |
| `NgsolveEddyBubbleHCurlBasis` | one-call NGSolve builder from high-order HCurl parent matrices/ports to an `EddyBubbleHCurlBasis` |
| `EddyTracePolynomialDim` | face trace mode count for simplex or tensor-product bridge enrichment |
| `EddyParentOrderLedger` | symbolic p-admissibility ledger: required parent order, p excess, and class-wise mode-count estimate |
| `NgsolveBridgeCycleCurrentBasis` | samples graph-cycle bridge currents as a VIM-compatible volume-current basis |
| `TopologyAwareHybridVIM` | production-facing tri-block assembly result: bulk EVRS, bridge-cycle current, and SIBC/surface-Omega blocks |
| `NgsolveTopologyAwareHybridVIM` | one-call NGSolve builder for topology classification, class-wise reduction planning, tri-block VIM assembly, and optional port RHS projection |
| `NgsolveEddyBubbleHybridVIM` | production station/stator builder: HCurl parent matrices/ports -> EVRS -> bulk/bridge/surface hybrid VIM |
| `AssembleHybridVIM` | assembles `R + sL + Z_s M_Gamma` on volume and surface bases |
| `HybridVIMSystem.block_rhs` | assembles full mixed right-hand sides from named volume/surface block projections |
| `HybridVIMSystem.diagnostics` | records the interaction backend name plus Hermitian/passive matrix checks for backend replacement |
| `HybridVIMSystem.schur_complement` | evaluates the IGTE mixed Galerkin block reduction `K_ss-K_sb K_bb^{-1}K_bs` |
| `ReducedInteractionMatrix` | accepts a pre-projected BEM/H-matrix interaction matrix |
| `ReducedPortAdmittance` / `ReducedPortImpedance` | evaluates the reduced external-port transfer function for p/n convergence and CLN fitting |
| `SIBCAdmittanceTail` | returns the leading `S sqrt(sigma/(mu s))` SIBC admittance tail |
| `SIBCSchurTerminationImpedance` | returns the Schur/Warburg scalar block `(s+d)/(K_SIBC sqrt(s))` |
| `SharedMeshMaterialModel` | carries the shared mesh, material labels, and `mu`/`nu`/`sigma`/SIBC coefficients seen by both branches |
| `NgsolveHDivMMMResponseReduction` | builds a protected physical-response basis plus energy-POD training enrichment for a supplied 3-D BDM or RT space |
| `NgsolveBDMHDivMMMResponseReduction` | constructs bare NGSolve `HDiv(order=p)` internally, locks the parent family as BDM, and builds the 3-D magnetic response reduction |
| `NgsolvePlanarHarmonicPorts` | builds 2-D real/imaginary harmonic-gradient training fields |
| `NgsolvePlanarHDivMMMResponseReduction` | applies the protected response-POD construction to `PlanarDemagBody`; BDM is the production default and `rt=True` is explicit |
| `MagnetizationCurrentCoupling` | builds the rectangular HDiv-magnetization / eddy-current coupling `int M_i dot B[J_j] dV` |
| `CoupleHDivMagnetizationToEVRS` | packages the rectangular HDiv-MMM / EVRS eddy-current coupling as a named mixed system |
| `CoupledHDivHybridVIMSystem` | HDiv-MMM coupled to a full hybrid HCurl-VIM eddy system: bulk EVRS, bridge-cycle, and surface-Omega/SIBC blocks |
| `MixedGalerkinHDivHybridVIMSystem` | frequency-specific exact Schur system obtained by eliminating only adjacency-class bulk eddy-bubble blocks from the complete HDiv/HCurl operator |
| `HCurlVIMHDivMMMSolution` | one-frequency mixed result with parent-HCurl `T`, sampled bulk/bridge/SIBC currents, sampled HDiv magnetization, Biot-Savart eddy `B`, loss, port response, and residual diagnostics |
| `CoupleHybridVIMWithHDivMMM` | builds the HDiv-MMM / hybrid HCurl-VIM coupling from a `HybridVIMSystem` and its current bases |
| `HCurlVIMHDivMMMSystem` | production-facing name for the mixed HCurl-VIM eddy branch and HDiv-MMM magnetic branch |
| `CoupleHCurlVIMWithHDivMMM` | production-facing constructor that records the optional shared mesh/material model |
| `CoupleEddyBubbleHCurlBasisWithHDivMMM` | direct HDiv-MMM coupling constructor for an `EddyBubbleHCurlBasis` |
| `NgsolveHCurlVIMHDivMMM` | one-call NGSolve production builder from HCurl parent matrices and an HDiv magnetization basis to the mixed MMM/VIM system |
| `NgsolveBDMEddyBubbleVIM` | complete production builder: shared mesh/material registry -> BDM response reduction + topology-aware eddy bubbling -> coupled frequency-domain VIM |

The important architectural choice is that the algebra is VIM.  The interaction
matrix is a replaceable VIM backend.  The same reduced basis can use:

```text
sampled dense Laplace VIM kernel
ngsolve.bem LaplaceSL projection as a VIM backend
Radia in-tree HACApK-compressed VIM projection
future HDiv-VIM material kernel coupling
```

Every production backend should pass the same `HybridVIMSystem.diagnostics`
gate: the reduced `R`, `L`, and `M_Gamma` blocks must be Hermitian and
positive semidefinite to numerical tolerance, and the backend name should be
recorded in validation JSON before CLN fitting or motor-coupled solves.

The port-level observable is

```text
Y_r(s) = B_r^* Z_r(s)^-1 B_r,
```

available through `ReducedPortAdmittance`.  This is the quantity to compare
under `p` refinement, Krylov/CLN depth changes, and backend replacement.

The current HDiv coexistence layer is intentionally rectangular.  The eddy
branch remains a current/impedance system, while the magnetic branch remains an
HDiv magnetization/demag system.  The shared block is the energy coupling

```text
C_ij = int_{Omega_m} M_i · B[J_j] dV.
```

This keeps the first implementation compatible with the existing HDiv-VIM
operator and avoids pretending that `HCurl` current DoFs and HDiv magnetization
DoFs are the same unknown.

The production-facing object for this is `CoupledHDivEVRSSystem`.  It records
the HDiv magnetization basis, the sampled EVRS current basis, the rectangular
coupling block, and optionally the reduced eddy-current impedance system.  This
is the intended bridge toward an HDiv-MMM magnetic branch coupled to an
eddy-bubble-eliminated HCurl eddy branch.

The production-facing eddy basis is now `EddyBubbleHCurlBasis`.  It is built
from an NGSolve parent space by `NgsolveEddyBubbleHCurlBasis`, which performs
the response-basis construction, samples `J=curl(T)`, and stores the
topology-aware eddy-bubbling decomposition.  The basis can then assemble its
own VIM block or couple to HDiv-MMM:

```python
eddy = vim.NgsolveEddyBubbleHCurlBasis(
    mesh, hcurl_fes, stiffness, mass, ports,
    steps=22,
    conductive_materials="cond",
    parent_order_ledger=vim.EddyParentOrderLedger(
        bulk_degree=4,
        bridge_trace_degree=0,
        surface_current_degree=2,
    ),
)
eddy_system = eddy.assemble_vim(sigma=5.8e7)
mixed = eddy.couple_hdiv_mmm(hdiv_magnetization_basis, eddy_system=eddy_system)
```

When the surface and bridge classes are retained, use the tri-block object
returned by `NgsolveTopologyAwareHybridVIM`.  It now has the same HDiv-MMM
coupling hook, but preserves the internal bulk/bridge/surface blocks:

```python
hybrid = vim.NgsolveTopologyAwareHybridVIM(
    mesh, hcurl_fes, response_vectors, surface_grad_modes,
    sigma=5.8e7,
    conductive_materials="cond",
)
mixed = hybrid.couple_hdiv_mmm(hdiv_magnetization_basis)
```

The resulting `CoupledHDivHybridVIMSystem` has one magnetic block and the full
hybrid eddy block.  The coupling matrix columns follow the eddy system order,
so Schur reduction can eliminate the eddy branch while retaining the distinction
between bulk EVRS, conductor-cycle bridge currents, and SIBC surface currents.

For production workflows that start from HCurl matrices rather than a
precomputed EVRS basis, use the one-call builder:

```python
multipole_ports = vim.NgsolveHDivRegularSolidHarmonicPorts(
    mesh,
    max_degree=3,
)
mixed = vim.NgsolveBDMEddyBubbleVIM(
    mesh,
    hcurl_fes,
    stiffness,
    mass,
    ports,
    surface_grad_modes,
    hdiv_order=1,
    mu_r=1001.0,
    external_fields=(H_stator_x, H_stator_y),
    external_names=("H_stator_x", "H_stator_y"),
    training_fields=multipole_ports,
    magnetic_materials="iron",
    steps=22,
    sigma=5.8e7,
    conductive_materials="cond",
    port_vector_potentials=port_vector_potentials,
)
solution = mixed.solve_frequency_eddy_bubbled(100.0)

T_parent = solution.parent_t_coefficients
J_bulk = solution.current_samples("bulk")
J_bridge = solution.current_samples("bridge")
K_sibc = solution.current_samples("sibc")
B_eddy = solution.eddy_flux_density(target_points)
M_hdiv = solution.sampled_magnetization
M_parent = solution.parent_magnetization_coefficients
```

The production builder records the neighboring-material roles explicitly as
`volume -> bulk`, `volume1 -> bridge`, and `surface -> sibc`.
`solve_frequency_eddy_bubbled()` obtains its keep/eliminate partition from
this role map: conductor-conductor cycle modes and conductor-air SIBC modes are
protected, while only the ordinary bulk candidate is eliminated.  It does not
apply one uniform basis rule to every element or face.

The mixed Galerkin transformation is formed from the complete coupled matrix,
not from the HCurl block in isolation.  For retained coordinates
`k = [HDiv, bridge, SIBC]` and eliminated bulk coordinates `b`, it constructs

```text
P_trial = [-A_bb^-1 A_bk; I]
P_test  = [-A_bb^-T A_kb^T; I]
P_test^T A P_trial = A_kk - A_kb A_bb^-1 A_bk.
```

Thus HDiv self interaction, both HDiv/HCurl coupling blocks, the projected RHS,
and the HCurl field lift are transformed together.  A nonzero eliminated RHS
uses the affine correction `x_b^(0) = A_bb^-1 f_b`, so reconstruction is an
exact solution of the original coupled reduced system rather than only a
homogeneous-basis projection.

`NgsolveBDMEddyBubbleVIM` first calls
`NgsolveBDMHDivMMMResponseReduction`, which creates bare NGSolve
`HDiv(order=hdiv_order)` and records `parent_family="BDM"`.  It then calls
`NgsolveEddyBubbleHybridVIM`: the eddy branch builds the response basis,
samples `curl(T)`, adds bridge-cycle and SIBC/surface-Omega blocks, assembles
the hybrid VIM impedance, and finally forms the HDiv-MMM coupling block.  The
lower-level `NgsolveHCurlVIMHDivMMM` remains available for explicit RT studies
or externally constructed magnetic bases.

The BDM-MMM parent demagnetizing operator remains the C++
`_ChargeGramHMatrix` HACApK path.  Mixed Galerkin acts after that parent
operator and the HCurl response basis have been projected.  The final coupled
Schur matrix is deliberately dense because its dimension is the retained HDiv
plus protected HCurl mode count, not the high-order parent DoF count.

For production naming the same concrete object is also exported as
`HCurlVIMHDivMMMSystem`, with constructor `CoupleHCurlVIMWithHDivMMM`.  This is
the name to use when describing the final method rather than the intermediate
EVRS research construction.

The mixed reduced system is written as

```text
[ A_m        alpha K ] [ m ] = [ f_m ]
[ beta K^*   Z_e(s) ] [ z ]   [ f_e ]
```

where `m` are HDiv-MMM magnetization coordinates, `z` are retained HCurl-VIM
eddy coordinates, `A_m` is the magnetic material/demag operator, `Z_e(s)` is
the reduced eddy impedance, and `K_ij = int M_i · B[J_j] dV`.  The API exposes
this algebra through `mixed_operator`, `solve`, `schur_magnetic_operator`,
`schur_eddy_operator`, and `port_admittance` methods on
`HCurlVIMHDivMMMSystem`.  The scale factors `alpha` and `beta` are explicit so
the final time/frequency convention can be selected without changing the
basis construction.

`NgsolveHDivMMMResponseReduction` first solves
`(M/(mu_r-1) + N) m_j = M H_j` for the physical and training fields with
Radia's native mass-Riesz C++ CG.  Physical stator/rotor response snapshots are
energy-orthonormalized and protected.  Regular-solid-harmonic or angle-sampled
training responses are projected into their magnetic-energy complement and
then POD-compressed.  Thus `max_modes` may discard enrichment, but it cannot
discard an operating excitation.

For regular solid harmonics through scalar degree `L`, the real port count is
`sum_(l=1)^L (2l+1) = L(L+2)`: 3, 8, and 15 ports for `L=1,2,3`.  Their H
fields have polynomial degree `L-1`, so the parent-space admissibility rule is
`HDiv order >= L-1` when exact polynomial representation is required.  The
response rank is selected from the saved energy-error curve, not from `L`
alone.

The lower-level `NgsolveHDivMMMReduction` applies the parent HDiv mass matrix and Radia
`DemagOperator` column by column, then forms `Q^* M Q`, `Q^* N Q`, and
`A_m = Q^* M Q / (mu_r - 1) + Q^* N Q` without constructing a dense parent
matrix.  Its `external_field_rhs` method projects new stator or rotor-angle
fields without rebuilding the demagnetizing operator.  `solve_frequency` uses
that stored magnetic operator and projected excitation RHS.  It derives the half-space skin
impedance from the stored conductivity unless an explicit SIBC impedance is
supplied.  The reduced complex solve uses the native row-major
`_HybridVIMSolve` kernel when `_radia_pybind` is available, with NumPy retained
as a source-level fallback.

### 6.1 End-to-End p=6 Production-Path Smoke

`validation_test/cln/hcurl_vim_hdiv_mmm_end_to_end.py` exercises the complete
path on a notched conductor with a p=6 HCurl parent, conductor-cycle bridges,
exterior-only SIBC modes, a response-adapted HDiv magnetization basis, and two
excitation ports.  The saved native-kernel smoke gives:

| parent HCurl DoF | Krylov | current-Gram bulk | bridge cycles | SIBC | total eddy modes |
|---:|---:|---:|---:|---:|---:|
| 3557 | 10 | 8 | 14 | 3 | 25 |

The production-default run uses NGSolve's bare `HDiv(order=1)`, hence BDM1,
with harmonic ports through degree two.  At 100 Hz and 10 kHz, the full mixed
RHS-scaled residuals are `4.69e-16` and `4.60e-16`, while the corresponding
backward errors are `2.58e-20` and `2.50e-18`.  It reconstructs all 3557 parent-HCurl
coordinates and all 267 parent-BDM1 coordinates; the magnetic branch reduces
267 DoFs to 8 modes.  The maximum snapshot residual is `9.91e-11` and the
full-rank training-response energy error is `3.54e-11`.  The same run evaluates
the summed bulk/bridge/SIBC eddy-current field at an exterior probe and records
`1.12e-4 T`, exercising the physical-field reconstruction rather than only the
reduced matrix solve.

The current-Gram quotient first removes 2 `curl(T)`-null directions and leaves
8 bulk modes whose Gram identity defect is `1.49e-14`.  The adjacency-role map
classifies `volume` as ordinary bulk, `volume1` as a
conductor-conductor cycle bridge, and `surface` as conductor-air SIBC.  The
mixed Galerkin step consequently eliminates 8 bulk coordinates but retains
all 8 BDM1, 14 bridge, and 3 SIBC coordinates, reducing the complete 33-mode
system to 25 modes.  Across 100 Hz and 10 kHz, the maximum full-coupled Schur
error is `1.11e-16`; direct-solve port response and Joule loss are reproduced
to `6.98e-16` and `8.14e-16`.  The parent HDiv demagnetizing action is supplied
by HACApK `_ChargeGramHMatrix`; only the small retained Schur solve is dense.

The companion `hcurl_vim_hdiv_mmm_bdm2_smoke.json` uses degree-three harmonic
ports and reduces 738 BDM2 DoFs to 15 modes.  Actual Raviart--Thomas is kept as
an explicit comparison, generated with `--hdiv-family rt`; RT1 and RT2 have
369 and 942 parent DoFs on this mesh.  Do not infer RT from the `HDiv` class
name.  Motor-specific training should replace or augment the generic harmonics
with actual slot and rotor-angle fields.

| HDiv parent | parent DoF -> modes | coupled modes -> retained | max Schur error | max port error | max loss error |
|---|---:|---:|---:|---:|---:|
| BDM1 | 267 -> 8 | 33 -> 25 | 1.110e-16 | 6.982e-16 | 8.140e-16 |
| BDM2 | 738 -> 15 | 40 -> 32 | 8.600e-17 | 4.238e-16 | 1.252e-15 |
| RT1 | 369 -> 8 | 33 -> 25 | 7.850e-17 | 9.090e-16 | 2.075e-15 |
| RT2 | 942 -> 15 | 40 -> 32 | 1.147e-16 | 6.657e-16 | 1.647e-15 |

The bulk block condition number was near `1e18` before the current-Gram
quotient.  The v8 complete coupled matrix is conditioned at `2.52e6` for 100 Hz
and `2.53e4` for 10 kHz.  Backward error, exact-Schur diagnostics,
reconstructed fields, port response, and Joule loss now agree simultaneously.

### 6.2 Planar BDM/RT Corner Smoke

`validation_test/cln/planar_hdiv_mmm_response_smoke.py` applies the same
construction to a 2-D L-shaped body.  The production-default BDM1 path reduces
42 DoFs to 4 modes and reproduces the parent Hx/Hy magnetization near the
re-entrant corner to `4.19e-12`.  BDM2 reduces 93 DoFs to 6 modes.  Explicit
RT1/RT2 comparisons have 62/123 parent DoFs.  These are response-span
reconstruction checks, not claims that a fixed rank covers every excitation.

### 6.3 Shared Mesh and Material Coefficients

The HCurl-VIM and HDiv-MMM branches may share the same NGSolve mesh and the same
material registry, but they should not share finite-element unknowns.  The
shared object is the region-labelled material model:

```text
mesh/material tags
  -> nu(x) or mu(x) for the HDiv-MMM magnetic branch
  -> sigma(x) and optional SIBC parameters for the HCurl-VIM eddy branch
```

The permeability is therefore shared by construction when both branches query
the same material registry on the same mesh.  It is not a scalar factor applied
after the eddy reduction.  The correct ordering is:

```text
build full material operators on mesh/regions
  -> apply T -> J = curl(T), HDiv magnetization, and surface-Omega maps
  -> project to retained EVRS/HCurl-VIM coordinates
  -> couple to HDiv-MMM through C_ij = int M_i · B[J_j] dV
```

This order preserves symmetry/passivity: the eddy bubble elimination sees the
same material operator that the full high-order problem would have seen.

The Python bridge records this contract with `SharedMeshMaterialModel`.  It is
intentionally light: `mu`, `nu`, `sigma`, and SIBC entries may be scalars,
per-region dictionaries, or backend coefficient functions.  Assembly code owns
the interpretation; the registry owns the fact that both branches are seeing
the same coefficients on the same mesh.

The dense complex Schur complement and the scalar SIBC/Warburg helper functions
are backed by Radia C++ kernels through `_radia_pybind` when the extension is
available, with NumPy fallbacks kept for source-level experimentation.

NGSolve's own `LOCAL_DOF` static condensation remains useful as an implementation
accelerator, but it is not the research concept called an eddy bubble here.  The
research reduction is EVRS/EID elimination: keep the response coordinates that
the external field can see, and eliminate the eddy-invisible directions.

---

## 7. What Must Be Proved or Demonstrated

For an IGTE-level contribution, the strong path is:

1. Define the high-order `HCurl` parent space and response-visible subspace.
2. Prove that the reduced VIM matrices preserve reciprocity and passivity.
3. Relate the block Krylov basis to CLN/Cauer moment matching.
4. Explain that `p` is parent-space resolution and `n` is circuit order.
5. Show that the surface-Omega/SIBC branch is another passive boundary branch.
6. Use small numerical examples only to illustrate p-vs-n convergence, VIM
   backend replacement, and the SIBC boundary branch.

The main message should be:

```text
High-order basis functions are not kept because the reduced model is large.
They are used so that the invisible DoFs can be eliminated after the correct
response space has been exposed.
```

---

## References

- A. Kameari et al., "Cauer Ladder Network Representation of Eddy-Current
  Fields for Model Order Reduction Using Finite-Element Method," IEEE
  Transactions on Magnetics, 2018.
- Y. Shindo et al., "Dynamical Model of an Electromagnet using Cauer Ladder
  Network Representation of Eddy-current Fields," IEEJ Journal of Industry
  Applications, 2018.  https://doi.org/10.1541/ieejjia.7.305
- O. Biro and N. Koster, "Generating a Cauer Ladder Network Representation of
  Eddy Current Fields Using Scalar Potentials," IEEE Transactions on Magnetics,
  2022.  https://doi.org/10.1109/TMAG.2022.3171079
- S. Yuferev and L. Di Rienzo, "Surface Impedance Boundary Conditions in Terms
  of Various Formalisms," IEEE Transactions on Magnetics, 2010.
- A. Bendali et al., high-order impedance-boundary/asymptotic treatments of
  magnetic skin effect and SIBC formulations.
