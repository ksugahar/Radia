"""Curated knowledge distilled from the IEEJ electromagnetic reports.

This is a decision-oriented synthesis of the technical-report archive covering
2005 through 2025.  It is deliberately not a raw bibliography dump: the
useful output is the method selection, limitation, and validation advice that
can guide Radia and NGSolve work.  The reports discuss both quasi-static and
wave problems; the current Radia policy remains MQS/Darwin with a Laplace
kernel, so wave/PML material is recorded as historical context only.
"""

OVERVIEW = r"""
# IEEJ electromagnetic technical reports: engineering map

The report series describes a consistent workflow for practical
electromagnetic CAE:

1. Choose the physical approximation and the observable before choosing a
   discretization.  In Radia, the default target is magneto-quasi-static (MQS)
   or Darwin behavior, not a full-wave Helmholtz solver.
2. Use NGSolve for finite-element plumbing: spaces, orientations, Piola maps,
   curved geometry, quadrature, weak forms, and grid-function evaluation.
3. Add a Radia-specific method only where it represents engineering physics or
   a reduced model that the general FEM stack does not provide.
4. Reduce repeated solves with model reduction, equivalent circuits, domain
   decomposition, or a surrogate, but re-evaluate the final design with the
   governing field model.
5. Treat force, loss, and optimization values as evidence packages.  Record
   the model, material law, operating point, mesh, solver, units, and error
   check together with the scalar result.

## Cross-cutting lessons

- Open-boundary accuracy is a formulation and boundary-condition problem, not
  merely a request for a larger air box.
- Material loss cannot be separated from the material scale and waveform.
  Classical eddy, anomalous eddy, hysteresis, stress dependence, lamination,
  and powder-core effects have different validity ranges.
- A reduced model is most useful when its ports, expansion points, states, and
  error range are explicit.  A fast model without a validity envelope is not a
  trustworthy surrogate.
- AI is a search and approximation layer around physics, not a replacement for
  the final field solve.  Active learning should spend new solves where the
  surrogate is uncertain or out of distribution.
- Force extraction needs a separate numerical error budget.  A field solution
  can look converged while a small force is still dominated by air-region or
  stress-integration error.
- High performance is application-dependent.  A preconditioner or Krylov
  method that wins for one mesh and coefficient distribution is not a universal
  default.

## Radia implementation posture

Use the existing NGSolve-native Kelvin, CLN/MOR, material, force, motor, and
HDiv/VIM knowledge tools for detailed API guidance.  This topic is the report
level map: it explains why a method is chosen and what must be checked before
promoting it into a panel, MCP workflow, or public result.
"""

OPEN_BOUNDARY = r"""
# Open-boundary methods: report comparison and current default

The early open-boundary reports compare four families.  Their historical
results are useful as benchmarks and as warnings about approximation order.

## Strategic Dual Image (SDI)

SDI combines boundary conditions with image fields so low-order multipoles
cancel at a finite artificial boundary.  The report examples give these
orders and choices:

- 2D circle: a 1:1 average of Dirichlet and Neumann conditions works for
  arbitrary multipoles in the reported Laplace setting.
- Axisymmetric ellipsoid, aspect ratio about 1.815: a 1:1 average is valid
  through the dipole term.
- Axisymmetric ellipsoid, aspect ratio about 1.390: weights about 0.424:0.576
  reach the reported quadrupole order.
- 3D sphere: the reported Neumann tangential-field and Dirichlet normal-field
  average is 1:2 and is valid through the dipole term.

These numbers are a method-selection reference, not a license to mix boundary
traces or field conventions.  Validate against an analytic multipole or a
larger-domain solve before using SDI in a production workflow.

## Kelvin transformation

Kelvin inversion maps an unbounded exterior to a bounded computational region.
The report examples use a 2D circle and 3D/axisymmetric spherical or ellipsoid
constructions.  In the current Radia implementation, use the NGSolve-native
Kelvin convention and its periodic boundary identification.  Do not replace a
failed Kelvin run by changing the material factor from memory: first check
matching boundary DOFs, orientation, the point representing infinity, gauge or
GND treatment, integration order, and refinement near the inversion center.

The report's historical axisymmetric extension uses a virtual permeability
that varies approximately with the inverse square radius.  It reaches the
reported dipole or quadrupole accuracy only for the stated geometry and
parameter choices.  The current 3D spherical Radia convention is the
canonical route; the report construction is a useful independent benchmark.

## IABC and shell approximations

IABC places one or more virtual-media shells around the physical region and
chooses shell parameters from continuity and the outer boundary condition.  In
the Laplace examples, n shells correct the multipole expansion through the
reported 2n-pole order.  The same idea was extended in the reports to
Helmholtz/high-frequency settings, but that is outside Radia's current
Laplace-only MQS policy and must not be silently imported into an MQS solver.

## Equivalent-theorem perturbation correction

Another report route solves a finite-domain FEM problem, reconstructs the
boundary field, and applies a first-order correction from an equivalent
theorem.  It can be much cheaper than a second FEM solve and is exact for the
single-turn example reported, but is only approximate in general.  Hybrid
FEM-DBCI variants iterate boundary data and therefore pay for multiple FEM
solves.

## Selection rule

1. For Radia production: use the maintained Kelvin/NGSolve path when the
   geometry and formulation satisfy its convention.
2. For a new open-boundary method: begin with a multipole benchmark and record
   the first omitted multipole, boundary geometry, and convergence trend.
3. Use SDI, IABC, or perturbation correction as a comparison lane unless a
   dedicated validation establishes their error envelope for the application.
4. Never infer full-wave radiation capability from the historical report
   chapters.  The present kernel and physical policy are quasi-static.
"""

MATERIALS = r"""
# Materials, hysteresis, and loss

The reports repeatedly warn that a material model is a coupled choice of
geometry scale, waveform, frequency, temperature/stress state, and observable.

## Laminations and anisotropic homogenization

Laminated cores can be represented as an anisotropic bulk material when the
lamination scale is below the field-resolution scale.  The effective model
must still distinguish the normal and in-plane directions and must be checked
against a resolved or analytic reference for the intended frequency range.
Do not treat a homogenized lamination as an isotropic BH curve by default.

## Stress-dependent magnetization

The 2023 report describes a multi-domain particle model with random easy axes
and energy contributions from Zeeman, crystal-anisotropy, magnetoelastic, and
demagnetizing terms.  A play-type pinning field can represent the history
dependence.  An equivalent stress that only uses a scalar magnitude can
underestimate off-axis loading; an energy-based stress invariant is more
physical.  A stress-coupled model therefore needs the stress tensor, magnetic
field direction, and history, not only one scalar stress.

## Anomalous eddy current

The grain-orientation model in the report captures low-induction anomalous
loss in a nonlinear A-phi setting, but at high induction it approaches the
classical eddy behavior and misses part of the anomalous contribution.  A
single fitted coefficient should not be extrapolated across that transition.

For the excess-eddy coefficient, the accompanying integral-equation paper
identifies its dependence on flux density directly from measured excess loss,
without imposing an arbitrary functional fit.  This is a good source for the
Radia field-separation and Bayesian identification lane: identify the
coefficient, carry its uncertainty, then re-evaluate the loss in the Cauer or
field model.

## Powder cores

The isolated-particle cell model in the report underpredicts eddy loss by
roughly a factor of five in the cited case.  Particle-size, particle-shape,
gap, domain-wall, and anomalous-loss variability matter.  Loss rises as the
microstructural variability is represented more realistically.  A powder-core
model must therefore expose its microstructural assumptions instead of hiding
them in one effective conductivity.

## Play plus Cauer

The report combines a play model for DC hysteresis with a Cauer ladder for
skin-effect and dynamic minor-loop behavior.  One ladder stage can overestimate
PWM minor-loop behavior; two stages improve the reported case.  A useful
workflow is to identify static hysteresis first, identify the dynamic ladder
from frequency-dependent loss, and validate an unseen waveform rather than
only the calibration loop.

## Darwin approximation

The Darwin model retains the Coulomb/longitudinal displacement-current effect
needed for parasitic capacitance while omitting the solenoidal/transverse wave
displacement current.  It is therefore an MQS/intermediate-frequency model,
not a radiation solver.  The A-Phi system requires a Coulomb-type gauge and
can become nonsymmetric if written naively; a symmetrized system and a scalar
gauge correction improve compatibility with symmetric solvers.  Keep the
solenoidal electric field in postprocessing even though its displacement
current is omitted from the source term.

## Material-model gate

Before a production run, state: material scale, BH/hysteresis law, loss
decomposition, lamination or particle geometry, stress/temperature, waveform,
frequency, and calibration range.  Then run at least one independent check:
analytic loss, resolved microstructure, an unseen waveform, or a measured
curve.  The result should expose which terms are modeled and which are not.
"""

MOR = r"""
# Model reduction, CLN, and circuit coupling

The reports use model reduction to turn repeated field solves into small
multiport networks.  The reduction is valuable only when ports and validity
conditions remain explicit.

## Multiport CLN

The multiport Cauer ladder represents a field region by port voltage/current
relations.  Multiple expansion points extend one reduced model over a wider
frequency range than one local expansion.  A useful record contains the port
ordering, expansion points, retained order, passivity constraints, source and
load convention, and a comparison against the unreduced field solve.

## Nonlinear CLN

For saturation, the report parameterizes the ladder basis by a state such as
the first ladder current and obtains fixed-state bases from static field
solutions.  The nonlinear workflow sweeps the state, solves the field problem
at those states, and interpolates or selects the corresponding reduced model.
Do not use a linear CLN outside the reluctivity range used to construct it.

## Lanczos relationship

The report explains CLN as a circuit interpretation of Lanczos reduction for
two Hermitian operators.  The tridiagonal and diagonal structures are useful
for checking implementation consistency, but they do not remove the need to
check symmetry, positive-real/passive behavior, port units, and residual error.

## Identification and sensitivities

An R/L ladder can be identified from FEM or measured impedance under
nonnegative/passive constraints.  Adjoint sensitivities are preferable to a
blind global search when a parameter gradient is available.  Complex-valued
derivatives need an explicit real/imaginary convention; otherwise a fit can
look numerically good while its gradient has the wrong sign.

## Motion and skew

For a moving conductor with position-dependent inductance L(p), the induced
voltage includes the directional derivative of L along the motion.  The
co-energy force relation is useful as a reduced-model check, but the position
derivative must be taken consistently with the current and port convention.

For a skewed induction motor, a slice model that solves every slice
independently misses axial current continuity.  A continuity-enforcing model
improves the physics but may constrain spatial-harmonic amplitudes.  A
per-slice-harmonic model is most faithful in the report, at the cost of more
states.  The choice is a modeling tradeoff, not just a mesh option.

## Darwin reduction

The reduced Darwin system can be constrained and nonsymmetric.  The report
uses nonsymmetric Arnoldi reduction and solves symmetric subblocks by
back-substitution, followed by a scalar correction enforcing the Coulomb-type
gauge.  Do not assume a symmetric CLN implementation applies unchanged to
this system.

## MOR validation gate

For every reduced model, compare frequency, transient, saturation, motion, or
skew behavior at points not used for construction.  Report the error in the
observable of interest, not only a matrix norm.  Reject extrapolation beyond
the expansion/state/port envelope or fall back to the full NGSolve solve.
"""

OPTIMIZATION = r"""
# Optimization and AI-assisted design

The reports treat optimization as a solve-budget problem.  The central rule
for Radia is: a surrogate proposes; the field model decides.

## Surrogate and Bayesian optimization

Gaussian-process Bayesian optimization uses a mean and uncertainty to balance
exploitation and exploration.  An acquisition rule that includes the
uncertainty term is useful only when the uncertainty is calibrated.  Active
learning should request new field solves near high uncertainty, sharp response
changes, or suspected out-of-distribution designs.

## Topology and parameterization

Density methods explore broadly but can leave gray or non-manufacturable
regions.  Level-set methods produce clearer interfaces but depend more on the
initial topology.  A practical two-stage route is density or parameter search
for broad exploration followed by level-set refinement.  A hybrid design can
keep a parameterized permanent magnet while giving barriers or air regions
topological freedom.

## Neural methods

The report surveys CNN/transfer learning, autoencoder latent spaces, SHAP or
Grad-CAM interpretation, and PINN-style approaches.  The practical warning is
that material-distribution generalization is uncertain.  Preserve the design
variables, constraints, operating points, training domain, and uncertainty or
OOD score with the prediction.

## Optimization gate

1. Define the physical objective and sign/unit convention.
2. Define feasibility constraints, manufacturability, stress, thermal, and
   operating-point limits.
3. Generate a DOE or initial field-solve set that covers the domain.
4. Fit the surrogate and test held-out designs.
5. Re-evaluate every selected design with the governing field model and retain
   the discrepancy.
6. Expand the training set where the discrepancy or uncertainty is large.

An AI result without the final field re-evaluation is a candidate, not a
validated electromagnetic design.
"""

HPC = r"""
# High-performance FEM and domain decomposition

The report series covers multigrid, FMM, domain decomposition, adaptive
meshing, high-order FEM, space-time parallelism, and large-scale iterative
solvers.  The implementation lesson is to measure the whole solve and choose
the method for the coefficient distribution and hardware.

## Static-field DDM

The static FEM system uses a curl-curl operator with a Coulomb constraint.
Subdomain Schur complements and a coarse correction control global modes.
Reported BDD and BDD-DIAG variants can reduce iterations dramatically in a
shaft example, but BDD-DIAG was not uniformly better than diagonal scaling in
that report.  Keep the comparison empirical.

## SVD local solves

For singular or nearly singular local subdomain systems, an SVD-based
Moore-Penrose inverse with a documented singular-value threshold can provide a
robust local solve.  It pays a memory and setup cost, then reuses the factor
for multiple right-hand sides.  The threshold and residual must be recorded;
silently truncating modes changes the physical problem.

## Very large eddy-current systems

The 130-billion-DOF example shows that Krylov choice matters: COCR converged
in the reported case while COCG did not reach convergence within the stated
iteration limit.  This is evidence for a benchmark lane, not a universal
solver ranking.  Store iteration history, residual definition, memory, and
hardware with any performance claim.

## Space-time parallel FEM

PinSTFEM-like methods combine spatial domain decomposition with time
parallelism.  Communication and block preconditioning become more difficult
as spatial partitions increase, so time/space partitioning is application
dependent.  Compare wall time, memory, communication, and physical error,
not only the number of cores.

## NGSolve boundary

NGSolve owns FE orientation, element transformations, Piola maps, quadrature,
and weak-form assembly.  Radia's C++ or Python acceleration must consume those
objects and preserve the caller-owned TaskManager convention.  A hand-written
high-order basis reconstruction is not a safe substitute for NGSolve's mapped
evaluation.
"""

FORCE = r"""
# Force and torque: numerical evidence rules

The force chapters distinguish microscopic and macroscopic field
interpretations and show why a visually converged field is not enough for a
small force.

## Choose the force observable

- Use Lorentz force for an imposed-current conductor when the target is a
  current region.
- Use co-energy/virtual work or a weighted Maxwell stress in air for magnetic
  material bodies.
- State whether the field is microscopic or macroscopic before selecting a
  stress tensor.  Deformation and magnetoelastic force require free-energy
  terms beyond a purely magnetic stress calculation.

## Air-region error

Nodal force extraction can leave a nonzero force on air due to discretization.
The report's correction is to subtract the air-model force, which costs a
second solve, or to refine the force-integration region.  Levitation and other
small-force problems need local refinement in the force direction and the
transverse directions; inspect the residual air force, not only the target
force.

## Source-separation correction

With two source groups, compute the full force and the force with each source
group alone.  Subtracting the source-only contributions from the full result
can cancel first-order source-independent discretization errors.  N source
groups require N+1 solves.  Apply it to the dominant source groups first and
record the source definition, current convention, field frame, units, and
sign convention.

## Radia validation gate

For a new force route, include one analytic or symmetry check, a mesh
refinement check, and an independent extraction when feasible.  Record whether
the result is full revolution, per-radian, per-length, or a symmetry-sector
quantity.  A force table without the geometry, source, unit, and convention
metadata is not a reusable validation artifact.
"""

MOTOR = r"""
# Motor modeling and optimization

The 2022 and 2025 reports connect reduced circuits, field solves, loss models,
thermal limits, and optimization into a virtual-motor workflow.

## Model fidelity choices

Reluctance networks and one-dimensional magnetic circuits are useful for rapid
screening.  Frozen-permeability field models, multiport CLN, Cauer winding
models, and semi-implicit plant circuits extend the operating range.  Choose
the model from the output required: average torque, torque ripple, transient
current, strand loss, iron loss, temperature, or fault behavior.

Skew and induction-machine models must state how slice coupling and spatial
harmonics are handled.  Treating every slice independently can break current
continuity.  Winding simplifications can also change strand eddy loss even
when average copper loss appears plausible.

## Loss and thermal coupling

The report covers PM eddy loss, AC copper loss, sleeve and harmonic mitigation,
stress in high-tensile steel, iron and stray loss, thermal contact, and cooling.
Thermal predictions need explicit contact and cooling assumptions.  A motor
efficiency map without those conditions is not portable between designs.

## Virtual motor artifact

Every reduced or surrogate motor result should retain: geometry/topology,
material laws, winding representation, operating points, speed/frequency,
thermal boundary conditions, solver/reduction settings, and comparison points
against a field model.  This makes an AI-generated motor design auditable and
reproducible.
"""

IMPLEMENTATION_MAP = r"""
# Implementation map: where to continue in Radia MCP

This report synthesis is an index into the maintained implementation knowledge.
Use the detailed topic tool named below after reading the report-level gate.

| Report lesson | Radia/NGSolve knowledge to consult next |
|---|---|
| Kelvin open boundary and periodic matching | ``kelvin_transformation`` and ``kelvin_identify_post_hoc`` |
| FEM/BEM and open-boundary comparison | ``ngsolve_usage`` and ``fem_bem_schur`` |
| CLN, Cauer, expansion points, MOR | ``cln_sibc_orthogonal``, ``cln_3d``, ``bem_cln``, ``cln_sphere_dd`` |
| Hysteresis, FORC, play, field separation | ``magnetic_materials_hysteresis`` and the hysteresis-ID topics |
| HDiv/VIM and NGSolve-native spaces | ``hdiv_vim`` and ``ngsolve_usage`` |
| Force error and independent checks | ``force_validation`` |
| Motor reduced workflows | ``airgap_motor_workflow`` and ``motor_usage`` |
| Panels and result-bearing workflows | panel design and notebook-workbench knowledge |
| Release and deployment evidence | ``release_workflow`` and the repository release gate |

When a report idea becomes production code, promote it through the repository's
validation ladder: fast regression test, numerical validation artifact,
result-bearing documentation, then a notebook panel only after the method is
validated.
"""

SOURCES = r"""
# Source coverage and provenance

This synthesis was prepared from the technical-report archive in the IEEJ
electromagnetic-analysis collection:

- 2005-12: high-speed and large-scale numerical computation technology.
- 2014-09: advanced use technology, including open boundaries, materials,
  force, coupling, and unresolved practical issues.
- 2016-12: advanced electromagnetic analysis technology, including MOR,
  materials, open boundaries, and applications.
- 2020-03: advanced accuracy technology, including homogenization, loss, MOR,
  high-performance FEM, and optimization.
- 2022-12 and 2023-05: advanced application technology, with dedicated
  chapters on materials, CLN/MOR, optimization, HPC, force, and applications.
- 2025-04: motor modeling, loss, thermal/structural constraints, and AI/ML
  optimization.
- 2013-2015 Radiation Science collections: historical neighboring work on
  open-boundary and equivalent-theorem approaches; mixed wave/radiation
  chapters are background, not Radia's MQS implementation target.
- Integral-equation paper on identifying the flux-density dependence of the
  excess-eddy field coefficient: direct support for the field-separation and
  loss-identification guidance.

Administrative announcement files, layout templates, and a motivational essay
were cataloged but do not add solver or modeling rules.  Two legacy Word files
were not machine-readable; their technical draft content is represented by the
corresponding PDF version in the same archive.
"""

TOPICS = {
    "overview": OVERVIEW,
    "open_boundary": OPEN_BOUNDARY,
    "open-boundary": OPEN_BOUNDARY,
    "kelvin": OPEN_BOUNDARY,
    "materials": MATERIALS,
    "material": MATERIALS,
    "hysteresis": MATERIALS,
    "loss": MATERIALS,
    "mor": MOR,
    "cln": MOR,
    "cauer": MOR,
    "optimization": OPTIMIZATION,
    "optimisation": OPTIMIZATION,
    "ai": OPTIMIZATION,
    "hpc": HPC,
    "ddm": HPC,
    "force": FORCE,
    "torque": FORCE,
    "motor": MOTOR,
    "motors": MOTOR,
    "implementation_map": IMPLEMENTATION_MAP,
    "source_map": SOURCES,
    "sources": SOURCES,
}


def get_technical_reports_documentation(topic: str = "overview") -> str:
    """Return curated technical-report knowledge for one topic or ``all``."""

    key = str(topic or "overview").strip().lower().replace(" ", "_")
    if key in {"all", "everything", "*"}:
        sections = [
            OVERVIEW,
            OPEN_BOUNDARY,
            MATERIALS,
            MOR,
            OPTIMIZATION,
            HPC,
            FORCE,
            MOTOR,
            IMPLEMENTATION_MAP,
            SOURCES,
        ]
        return "\n\n".join(sections)
    if key in TOPICS:
        return TOPICS[key]
    available = ", ".join(sorted(TOPICS))
    return f"Unknown topic {topic!r}. Available topics: {available}, all"
