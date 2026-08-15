# EarlyTimes C++ API Design

**Date:** 2026-08-15

**Status:** Native tracking, transfer/Lie kernels, and HCurl GridFunction vector-potential input implemented

**Scope:** Target contract for migrating the useful EarlyTimes concepts into
the Radia C++ core, with pybind11 and standalone MEX bindings over the same
implementation.

## Implementation Status

Three implementation slices are now available in `radia::beam`:

1. A canonical six-dimensional, piecewise-constant variational kernel
   propagates `R`, `T`, and `U` through third order, preserves the public
   factorial convention, and attributes quadratic terms, direct cubic terms,
   local cubic cascades, and ordered upstream-to-downstream region-pair
   cascades.
2. A shared NGSolve adapter accepts a live real three-component vector
   potential `A` in `HCurl(order=p)`. NGSolve evaluates its native `curl(A)`
   at the center and eight transverse ring points per supplied reference
   station; Radia never projects it to HDiv and never constructs a regular-grid
   field map. The adapter fits source-free normal/skew multipoles, reports the
   retained FESpace class and order `p`, and feeds the native linear or
   higher-order map. Direct magnetic-flux-density GridFunctions remain an
   explicit `magnetic_flux_density` compatibility mode, not the default.
3. A dependency-free native tracking foundation provides validated SI
   `ParticleSpecies`, `ReferenceParticle`, and `CartesianState` values;
   inspectable zero and uniform electromagnetic fields; a relativistic
   Lorentz right-hand side with time, path length, or azimuth as the explicit
   independent variable; individually callable classical RK4 and relativistic
   Boris steps; and a fixed-step tracker returning every state, step, field
   sample, and invariant report. The Python and standalone MEX surfaces call
   this same C++ implementation. Pure-magnetic Boris tracking preserves
   momentum magnitude to roundoff in the analytic uniform-field gate.

Thin pybind11 and standalone MEX boundaries call these same C++ sources of
truth through `radia.beam.propagate_grid_function_linear_map` and
`radia.beam.propagateGridFunctionLinearMap`, respectively. Native field
sampling, `LorentzEquation.rhs`, one-step integration, and trajectory building
are likewise available independently in both languages. The compatibility
linear entry returns only `R`; the multipole entry fits through the requested
field order and propagates the declared nonlinear map.

General lattice composition, event handling, adaptive/Lie particle stepping,
real CAD design-parameter derivatives including edge boundary terms, and
ray-based map-validity checks remain later phases. No current API claims that
a local multipole profile alone identifies an edge-angle recommendation.

## Decision

The public namespace is `radia::beam`, not `radia::early_times`.
There is no EarlyTimes file-compatibility namespace. EarlyTimes is the
algorithm and validation corpus for the migration, not a public input/output
contract. Its global state, fixed arrays, cylindrical state layout, namelists,
and file-driven execution model are retired.

The numerical source of truth is C++. Python and MATLAB do not reimplement the
equations, integrators, closed-orbit search, transfer-map fitting, or optics.
pybind11 exposes inspectable C++ objects. MEX exposes the same components through
typed handles and stateless diagnostic commands.

The design deliberately does not provide only a monolithic `run()` entry point.
A user can inspect and call field evaluation, the equation right-hand side, one
integration step, an event function, a closed-orbit residual, a transfer-map
stencil, and each optics diagnostic independently.

The primary new deliverable is not another way to read `.sol` or to draw a
tracked orbit. Radia already has field serialization and point-evaluation
paths. The primary deliverable is a distributed optics analyzer that answers:

- where the solved three-dimensional field acts as a normal or skew
  quadrupole;
- how the local linearized dynamics accumulate into the endpoint transfer
  matrix;
- which body, edge, fringe, longitudinal-field, or alignment region creates a
  selected coupling term; and
- whether changing a real design parameter such as an edge angle strengthens
  or weakens the desired endpoint optics.

## Primary Deliverable: Distributed Transfer Attribution

For a reference trajectory `z0(s)`, the analyzer forms the full canonical
linearization

```text
A(s) = d f(s, z) / d z | z=z0(s)
d R(s, s0) / ds = A(s) R(s, s0),   R(s0, s0) = I.
```

It stores both `R(s, s0)` and `R(sf, s)`. The field is sampled directly from
the NGSolve `GridFunction` in the trajectory's local transverse frame. A local
least-squares multipole fit is the robust default for HDiv/HCurl fields because
raw elementwise derivatives may be discontinuous. The fit reports normal and
skew coefficients, sample radius, rank, conditioning, and residual. It never
constructs an intermediate regular-grid field map.

Two related results must remain distinct:

1. `k1_normal(s)` and `k1_skew(s)` identify where the local field has
   quadrupole character.
2. A transfer sensitivity density identifies where changing that character
   matters to the endpoint map.

For a named local generator component `A_c(s)`, its infinitesimal local-strength
sensitivity is

```text
D_c(s) = R(sf, s) A_c(s) R(s, s0).
```

For a real design parameter `p`, such as an entrance-edge angle, the required
quantity is instead

```text
dR(sf, s0)/dp = integral R(sf, s) (dA(s)/dp) R(s, s0) ds
                + boundary and reference-orbit terms.
```

`D_c(s)` is a Frechet sensitivity density, not a claim that the finite final
matrix has a unique additive decomposition. Reports therefore distinguish
local field role, local-strength leverage, and an actual design-parameter
derivative. This distinction is what turns a plot into a design decision such
as "increase the entrance edge angle" rather than merely noting a large field
gradient.

The default moving frame is a right-handed Bishop rotation-minimizing frame
seeded by an explicit design normal. Its discrete implementation is the
fourth-order double-reflection method of Wang, Juttler, Zheng, and Liu,
*ACM TOG* 27(1), 2008, DOI 10.1145/1330511.1330513, rather than repeated
normal-plane projection. Frenet framing is available as an explicit option.
The selected convention is stored in every result so frame rotation is not
misidentified as physical skew quadrupole or coupling.

For stations sampling a complete closed trajectory, the explicit periodic
mode converts the open-path RMF to the periodic constant-twist minimal-twist
frame by distributing the measured one-turn holonomy in chord arc length.
This follows the boundary-frame distinction of Farouki and Moon (2018), DOI
10.1007/s10444-018-9599-3.  It is not enabled implicitly: a symmetry-reduced
cell needs its known spatial symmetry in the closure map rather than a direct
last-to-first chord.

The existing planar machinery in `radia.isochronous_topopt` and
`radia.accelerator_magnet_topopt` is the first reusable specialization: it
already propagates piecewise-constant curvature and normal gradient into a
6-by-6 map and chains analytic matrix-exponential Frechet derivatives to field
responses. The new analyzer generalizes that chain to a solved 3-D field,
normal and skew quadrupoles, solenoidal terms, fringe/edge regions, arbitrary
reference trajectories, and station/interval-resolved output. The planar
combined-function result remains an independent regression oracle.

## Nonlinear Map And Attribution

The design decision for an edge or end pack must not be made from `R` alone.
Strengthening an edge may improve first-order focusing while worsening
geometric aberration, chromatic coupling, or amplitude-dependent behavior. The
canonical deviation `u` is therefore expanded through third order by default:

```text
f(s, z0 + u) = f0 + A(s) u
               + (1/2) F2(s)[u,u]
               + (1/6) F3(s)[u,u,u] + O(u^4)

u(sf) = R u0 + (1/2) T[u0,u0] + (1/6) U[u0,u0,u0] + O(u0^4).
```

The coefficient storage is symmetric in the input indices and the factorial
convention above is part of the public ABI. The variational tensors satisfy

```text
R' = A R
T' = A T + F2[R,R]
U' = A U + 3 F2[R,T] + F3[R,R,R],
```

where all repeated polynomial arguments are symmetrized. Production map
construction integrates these equations or an equivalent checked canonical
jet/TPSA representation. Finite-difference ray fitting remains an independent
regression route, not the production derivative engine.

The second-order endpoint tensor has the exact position density

```text
K_T(s) = R(sf,s) F2(s)[R(s,s0), R(s,s0)].
```

The third-order density contains both a direct local term and a cascade term:

```text
K_U_direct(s)  = R(sf,s) F3(s)[R,R,R]
K_U_cascade(s) = 3 R(sf,s) F2(s)[R,T].
```

The cascade cannot honestly be assigned to only one location because `T(s)`
already contains upstream nonlinearities. In addition to a one-dimensional
station summary, the analyzer therefore returns an ordered region-pair table
for dominant `upstream source -> downstream interaction` contributions. This
exposes effects such as two separated sextupole/fringe regions producing a
third-order coupling term.

Multipoles are reported in two origins:

- the declared magnet/design center, which describes manufactured field
  quality; and
- the actual reference orbit, which describes beam-effective multipoles and
  includes offset feed-down.

The report includes normal/skew sextupole and octupole profiles, longitudinal
derivatives, geometric and chromatic `T`/`U` entries, amplitude detuning and
resonance-driving diagnostics derived from an inspectable Lie/normal-form
factorization, and the raw tensors from which each diagnostic was computed.
No normal-form scalar replaces the raw map.

Every nonlinear report declares an amplitude and momentum validity domain. A
deterministic ray set over that domain is directly tracked through the original
field and compared with the truncated map. Recommendations fail closed when
the requested operating domain exceeds the measured map-error tolerance.

An edge-angle recommendation uses a joint objective such as

```text
J = ||W_R (R - R_target)||^2
    + ||W_T T||^2 + ||W_U U||^2
    + selected resonance/chromatic penalties.
```

The report returns `dJ/dp` and separate linear, quadratic, cubic, chromatic,
and resonance contributions. It may recommend increasing an edge angle only
when the first-order gain survives the configured nonlinear limits over all
required momenta and amplitudes.

## NGSolve Design Principles Applied

EarlyTimes should follow the useful parts of the NGSolve object model:

1. **Composable fields.** A magnetic field behaves like a vector-valued
   `CoefficientFunction`: it can be evaluated, composed, transformed, scaled,
   differentiated when supported, and inspected as an expression tree.
2. **Mesh and field stay distinct.** The production Lie/A-RK input is a real
   NGSolve `GridFunction` in `HCurl(order=p)`, representing the vector
   potential `A`.  The independent validation input is preferably the original
   HDiv-MMM magnetic-flux-density `CoefficientFunction`; an `HDiv`
   GridFunction is an optional projection check.  NGSolve owns point location
   and field evaluation; Radia never reconstructs finite-element basis data.
3. **Operators are explicit.** The equation of motion, one-step operator,
   one-turn map, closed-orbit residual, and transfer map are separate objects.
4. **State is visible.** Particle state, reference frames, accepted steps,
   rejected steps, field samples, events, residuals, and fit diagnostics are
   ordinary value objects, not hidden global arrays.
5. **Composition precedes orchestration.** `Tracker::Track` is a convenience
   composition of public lower-level operations, not the only way to execute
   them.
6. **No finite-element plumbing is duplicated.** NGSolve owns mapped points,
   curved geometry, element search, `CoefficientFunction`, and `GridFunction`
   evaluation. Radia owns beam dynamics and field-to-beam attribution.

NGSolve requires a mapped integration point for mesh-dependent fields.
Therefore `NGSolveGridFunctionField` retains the native HCurl `GridFunction`,
obtains the mesh through its `FESpace`, and reports a clear outside-mesh status
instead of silently extrapolating. Its default magnetic coefficient is the
NGSolve `GridFunctionCoefficientFunction::Deriv()` result. The input order is
read from the FESpace; the beam API does not accept a second, possibly
inconsistent `p` argument.

## Layering

```text
Python radia.beam                 MATLAB +radia/+beam
        |                                |
        +----------- pybind11 -----------+-- standalone radia_mex commands
                         |
                    radia::beam
       +-----------------+-------------------------------+
       |                 |               |               |
     Field            Dynamics        Tracking       Diagnostics
  composition         equation        steppers        and optics
       |                 |               |               |
       +-----------------+---------------+---------------+
                         |
              Radia field | NGSolve CF/GF | analytic FFAG
```

The core has no dependency on Python or MATLAB. The NGSolve adapter is a thin
native layer above the generic field contract. Configuration enters through
typed value objects, and result serialization remains outside the numerical
core.

## Units And Coordinates

All core values use SI units.

| Quantity | Unit |
|---|---|
| position, path length | m |
| time | s |
| kinetic momentum | kg m/s |
| electric field | V/m |
| magnetic vector potential | T m |
| magnetic flux density | T |
| charge | C |
| rest mass | kg |
| angle | rad |
| kinetic energy at API boundary | eV or J, named explicitly |

The production A input is a three-component real NGSolve `HCurl(order=p)`
`GridFunction` whose value is the vector potential in tesla-metres.  The Lie
and canonical A-RK routes evaluate A itself; they do not replace it with
`curl(A)`.  HDiv-MMM generates continuous A and B source
`CoefficientFunction` objects. A is conformingly projected to HCurl p=5.  The
independent Cartesian B-RK route may evaluate the original B CoefficientFunction
directly, while HDiv p=4 remains an optional projected-field comparison.  The
field spaces, mesh transformations, orientation, and evaluation remain owned
by NGSolve. The beam API does not build or consume a regular-grid field map.

Lie point evaluation locates the containing element with the HCurl
GridFunction's mesh and evaluates the vector potential itself at the mapped
point.  `track_hcurl_vector_potential_canonical_s` does the same for `A_s,A_y`
and recovers transverse derivatives from fourth-order centred differences of
those NGSolve-owned point values in the exact unexpanded Hamiltonian. It does
not ask HCurl for a nonconforming element-interior `Grad(A)` field.
`track_b_coefficient_cartesian_s` independently evaluates
HDiv-MMM's supplied B CoefficientFunction in Cartesian Lorentz RK, or an
HDiv GridFunction when explicitly testing that projection.  `curl(A)` is a
diagnostic consistency quantity, never the B-map validation substitute.
`compare_hcurl_lie_map_to_direct_rk` records Lie-versus-A truncation,
A-versus-B field-route discrepancy, and total Lie-versus-B discrepancy as
separate arrays.  The B route consults the A gauge only when converting
canonical to mechanical momentum at the entrance/exit.

Measured median-plane magnetic fields enter one stage earlier.  They define
`MeasuredMedianPlaneFieldTarget` rows for HDiv-MMM pole-topology optimization
at the physical probe locations.  Every accepted pole change is followed by a
complete three-dimensional field re-solve.  No B-spline/polynomial continuation
of a measured plane is accepted as off-plane B, HCurl A, Lie input, or RK input.
Only the A and B fields generated from the accepted physical magnet solution
cross the EarlyTimes boundary; held-out measurements may independently check
that solution.

For a fourth-order A-map, `sample_transverse_vector_potential` evaluates the
HCurl GridFunction on the full upper/lower local `(x,y)` patch.
The public `fourth_order_lie_map_from_hcurl_transverse` boundary certifies
HCurl order, `A_x=0`, design-orbit `A_s=A_y=0`, and the declared normal/skew
symmetry before its private degree-five jet recovery.  Arbitrary polynomial
coefficients are not a field-map input. Separate left/right patch fits expose
a derivative mismatch at the central `x=0` face without defining `A_y` by a
two-trace average. The Lie engine differentiates the complete six-dimensional
`R/T/U/V` and `f3/f4/f5` tensors in forward mode. Under normal median-plane
symmetry the first-order horizontal/vertical blocks vanish, while all
symmetry-allowed second- through fourth-order xy and momentum-offset cross
terms remain available for zero or nonzero optimization targets;
`differentiate_hcurl_transverse_lie_map` chains any supplied HCurl-DOF or
topology response through the complete fit-and-map calculation.

The map builder requires the Hamiltonian linear term to vanish within a
declared tolerance.  It returns the derivative of that term as the fixed-orbit
constraint Jacobian.  Optimization must either remain in this constraint's
null space or recompute and differentiate the design orbit; silently dropping
the affine orbit displacement is forbidden.

Two state types are public and must not be conflated:

```cpp
struct CartesianState {
    Vec3 position_m;
    Vec3 kinetic_momentum_kg_m_s;
    double time_s;
    double path_length_m;
};

struct CanonicalState6 {
    double x_m;
    double px_over_p0;
    double y_m;
    double py_over_p0;
    double sigma_m;
    double delta;
};
```

`CanonicalState6` is defined in a documented right-handed Frenet or parallel-
transported frame. The old ET `(r, theta, z, p_r, p_theta, p_z)` array is not a
public input format. Slope `(x, x', y, y')` values may be converted explicitly
for reporting, but are not the internal state.
Transfer maps use canonical coordinates so their symplectic defect is
meaningful.

## Core Value Types

```cpp
namespace radia::beam {

struct Vec3 { double x, y, z; };
struct Mat3 { std::array<double, 9> row_major; };

struct ParticleSpecies {
    double charge_c;
    double rest_mass_kg;
    std::string name;
};

struct ReferenceParticle {
    ParticleSpecies species;
    double kinetic_energy_j;
    double momentum_kg_m_s;
    double magnetic_rigidity_t_m;

    static ReferenceParticle FromKineticEnergyEV(
        ParticleSpecies species, double kinetic_energy_ev);
};

struct Frame3 {
    Vec3 origin_m;
    Vec3 ex;
    Vec3 ey;
    Vec3 tangent;
};

enum class DomainStatus { inside, outside, boundary, invalid };

struct FieldRequest {
    bool electric = false;
    bool magnetic = true;
    bool spatial_jacobian = false;
    bool source_breakdown = false;
};

struct FieldContribution {
    std::string source_name;
    Vec3 electric_v_m;
    Vec3 magnetic_t;
};

struct FieldSample {
    Vec3 electric_v_m;
    Vec3 magnetic_t;
    std::optional<Mat3> dB_dx_t_m;
    DomainStatus domain_status;
    std::vector<FieldContribution> contributions;
};

}  // namespace radia::beam
```

Every public struct validates finite values and physical constraints at its
constructor or factory boundary. No API interprets an unnamed six-element
array.

## Field API

```cpp
class Field {
public:
    virtual ~Field() = default;
    virtual FieldSample Evaluate(
        const Vec3& position_m,
        double time_s,
        const FieldRequest& request = {}) const = 0;
    virtual void EvaluateBatch(
        ArrayView<const Vec3> positions_m,
        double time_s,
        const FieldRequest& request,
        ArrayView<FieldSample> output) const;
    virtual FieldDescription Describe() const = 0;
    virtual std::vector<std::shared_ptr<const Field>> Inputs() const = 0;
};
```

For the default `NGSolveGridFunctionField`, `Evaluate` is conceptually:

```cpp
auto mesh = grid_function->GetMeshAccess();
auto magnetic_flux_density = grid_function->Deriv(); // native HCurl curl(A)
auto element = mesh->FindElementOfPoint(position, integration_point, true);
auto& transformation = mesh->GetTrafo(element, local_heap);
MappedIntegrationPoint<3, 3> mapped(integration_point, transformation);
magnetic_flux_density->Evaluate(mapped, result);
```

An invalid element search returns `DomainStatus::outside`. No nearest-cell,
clamping, or sampled-map fallback is allowed.

Required concrete fields:

- `ZeroField`
- `UniformField`
- `IdealRadialSectorFFAGField`
- `NGSolveGridFunctionField` retaining a three-component native
  `ngcomp::GridFunction` in `HCurl(order=p)` and its native curl coefficient;
  its space and mesh are retained through shared ownership. An explicit
  direct-B mode supports independent HDiv/VectorH1 comparisons.
- `SumField`, `ScaledField`, `TransformedField`, and `RegionField`
- `ErrorField` for alignment, strength, and measured-error terms

Composition remains inspectable:

```cpp
auto field = SumField({
    std::make_shared<NGSolveGridFunctionField>(solved_vector_potential),
    std::make_shared<TransformedField>(measured_error, alignment),
});

auto tree = field->Describe();
auto sample = field->Evaluate(point, time, {.source_breakdown = true});
```

`FieldDescription` is a serializable expression tree containing type, name,
units, coordinate frame, children, interpolation policy, and provenance. It is
available identically from Python and MATLAB.

## Lattice API

The lattice is data, not an implicit control flow encoded in element names.

```cpp
struct Aperture {
    virtual ApertureResult Evaluate(const CartesianState&) const = 0;
};

struct ThinOperator {
    virtual CanonicalState6 Apply(
        const CanonicalState6&, const ReferenceParticle&) const = 0;
    virtual LinearMap6 Linearize(
        const CanonicalState6&, const ReferenceParticle&) const = 0;
};

struct Element {
    std::string name;
    Transform3 placement;
    double design_length_m;
    std::shared_ptr<const Field> field;
    std::shared_ptr<const Aperture> aperture;
    std::vector<std::shared_ptr<const ThinOperator>> entrance_operators;
    std::vector<std::shared_ptr<const ThinOperator>> exit_operators;
    StepPolicy step_policy;
};

class Lattice {
public:
    const std::vector<Element>& Elements() const;
    FieldSample EvaluateField(const CartesianState&, const FieldRequest&) const;
    ElementLocation Locate(const CartesianState&) const;
    LatticeDescription Describe() const;
};
```

Drifts, magnets, septa, kickers, markers, hard edges, and fringe maps are
represented explicitly. A hard-edge kick is a `ThinOperator`; it is never a
side effect hidden inside field evaluation. `Lattice::Describe()` exposes the
ordered element list, placements, field expression trees, apertures, and step
policies.

## Dynamics And Integrators

```cpp
enum class IndependentVariable { time, path_length, azimuth };

class Equation {
public:
    virtual StateDerivative RHS(
        double independent_value,
        const CartesianState& state,
        EvaluationRecord* record = nullptr) const = 0;
    virtual InvariantReport Invariants(const CartesianState&) const = 0;
};

class LorentzEquation final : public Equation {
public:
    LorentzEquation(ReferenceParticle, std::shared_ptr<const Field>,
                    IndependentVariable);
};

class Stepper {
public:
    virtual StepResult Step(
        const Equation&, double independent_value,
        const CartesianState&, double step) const = 0;
    virtual StepperDescription Describe() const = 0;
};
```

Required steppers are:

- `ClassicalRK4`, matching the canonical EarlyTimes `rkclassic` route
- `LieDriftKick2`, matching the canonical second-order Lie route
- `Boris2`, the modern magnetic-field default because it preserves momentum
  magnitude for a pure magnetic field
- an explicit adaptive step controller that wraps a named embedded method and
  returns accepted and rejected step records

The public `RHS` and `Step` methods are mandatory. This is the primary guard
against black-box behavior: users can inspect the field and derivative at one
state, then reproduce every accepted trajectory step themselves.

## Tracking API

```cpp
struct TrackPlan {
    double start;
    double stop;
    StepSchedule steps;
    SampleSchedule samples;
    std::vector<std::shared_ptr<const Event>> events;
    TraceLevel trace_level = TraceLevel::accepted_steps;
};

struct StepRecord {
    double independent_value;
    double attempted_step;
    double accepted_step;
    bool accepted;
    CartesianState state_before;
    CartesianState state_after;
    std::optional<FieldSample> field;
    std::optional<StateDerivative> rhs;
    InvariantReport invariants;
    std::string element_name;
};

class Trajectory {
public:
    const std::vector<CartesianState>& Samples() const;
    const std::vector<StepRecord>& Steps() const;
    const std::vector<EventRecord>& Events() const;
    TrajectorySummary Summary() const;
};

class Tracker {
public:
    Tracker(std::shared_ptr<const Equation>, std::shared_ptr<const Stepper>);
    StepResult Step(double independent_value, const CartesianState&,
                    double step) const;
    Trajectory Track(const CartesianState&, const TrackPlan&) const;
};
```

Trace levels allow performance control, but summaries never suppress failure,
event, domain, or invariant information. The core returns data and does not
open files. JSON, MAT, and GMSH writers are separate exporters.

## Reference Orbit And Closed Orbit

Closed-orbit solving is decomposed into a map and a nonlinear solver:

```cpp
class OneTurnMap {
public:
    CanonicalState6 Apply(const CanonicalState6&) const;
    LinearMap6 Linearize(const CanonicalState6&, const DifferenceStencil&) const;
};

struct ClosedOrbitIteration {
    int iteration;
    CanonicalState6 state;
    CanonicalState6 residual;
    LinearMap6 jacobian;
    double residual_norm;
    double accepted_damping;
};

class ClosedOrbitProblem {
public:
    CanonicalState6 Residual(const CanonicalState6&) const;
    LinearMap6 Jacobian(const CanonicalState6&, const DifferenceStencil&) const;
    ClosedOrbitResult Solve(const CanonicalState6&, const NewtonOptions&) const;
};
```

`Solve()` returns every iteration. Random restart, scan, and correction are
separate strategies operating on the same public residual. They do not live as
flags in one global control array.

## Transfer Maps, Optics, And Coupling Attribution

```cpp
struct TransferMap {
    CanonicalState6 reference_input;
    CanonicalState6 reference_output;
    LinearMap6 first_order;
    Tensor6x6x6 second_order;
    std::optional<ThirdOrderMap6> third_order;
    CoordinateScales scales;
    FitReport fit;

    CanonicalState6 Apply(const CanonicalState6&) const;
    double SymplecticDefect() const;
};

class TransferMapBuilder {
public:
    DifferenceStencil Stencil(const CanonicalState6&) const;
    std::vector<TrackedRay> TrackStencil(const CanonicalState6&) const;
    TransferMap Fit(const std::vector<TrackedRay>&) const;
    TransferMap Build(const CanonicalState6&) const;
};

struct OpticsReport {
    LinearMap6 one_turn_map;
    CoupledTunes tunes;
    TwissReport twiss;
    DispersionReport dispersion;
    ChromaticReport chromatic;
    ResonanceReport resonance;
    double symplectic_defect;
};
```

`TransferMapBuilder::Stencil`, `TrackStencil`, and `Fit` are public so a map is
not an unexplained matrix. The fit report contains perturbation sizes, ray
states, rank, condition estimate, residuals, and rejected samples.

Ray-stencil fitting is the independent nonlinear-map route. The primary
first-order diagnostic uses the variational equation directly so it can retain
the position at which each effect enters the map.

Field-to-beam attribution is a Radia extension beyond the legacy ET output:

```cpp
enum class GeneratorComponent {
    drift,
    curvature,
    normal_quadrupole,
    skew_quadrupole,
    longitudinal_field,
    fringe_and_edge,
    electric,
    unclassified
};

struct GeneratorDecomposition {
    std::map<GeneratorComponent, LinearMap6> parts_per_m;
    LinearMap6 reconstruction_residual_per_m;
};

struct LinearizationStation {
    double path_length_m;
    Frame3 frame;
    MultipoleSpectrum multipoles;
    MultipoleFitReport multipole_fit;
    LinearMap6 generator_a_per_m;
    GeneratorDecomposition generator_parts;
    LinearMap6 r_from_start;
    LinearMap6 r_to_end;
    std::optional<Tensor6x6x6> cumulative_t;
    CouplingMetrics coupling;
    std::vector<FieldContribution> field_sources;
};

struct TransferSensitivityDensity {
    std::string parameter_name;
    std::vector<LinearMap6> density_per_m;
    LinearMap6 integrated_dR_dp;
    std::optional<double> objective_derivative;
};

struct IntervalSensitivity {
    double s_begin_m;
    double s_end_m;
    std::string element_name;
    LinearMap6 integrated_dR_dp;
    double signed_objective_derivative;
    double absolute_objective_leverage;
};

struct TransferAttributionReport {
    TransferMap endpoint_map;
    std::vector<LinearizationStation> stations;
    std::vector<IntervalSensitivity> ranked_intervals;
    LinearizationDiagnostics diagnostics;
};

struct NonlinearMapReport {
    LinearMap6 r;
    Tensor6x6x6 t;
    ThirdOrderMap6 u;
    std::vector<SecondOrderDensity> t_density;
    std::vector<ThirdOrderDensity> u_direct_density;
    std::vector<OrderedRegionPairContribution> u_cascade_pairs;
    NormalFormReport normal_form;
    MapValidityReport validity;
};

class TransferAnalyzer {
public:
    LinearizationStation LocalLinearization(
        const ReferenceSample&, const LinearizationOptions&) const;
    TransferAttributionReport Analyze(
        const Trajectory& reference,
        const AttributionOptions&) const;
    TransferSensitivityDensity LocalStrengthSensitivity(
        const TransferAttributionReport&,
        GeneratorComponent,
        const MatrixObjective&) const;
    TransferSensitivityDensity DesignSensitivity(
        const TransferAttributionReport&,
        const DesignParameterDerivative&,
        const MatrixObjective&) const;
    NonlinearMapReport AnalyzeNonlinear(
        const Trajectory& reference,
        const NonlinearAttributionOptions&) const;
    NonlinearDesignSensitivity NonlinearSensitivity(
        const NonlinearMapReport&,
        const DesignParameterDerivative&,
        const NonlinearMapObjective&) const;
};
```

`GeneratorDecomposition` must reconstruct `generator_a_per_m` within a reported
tolerance; terms which cannot be assigned without convention go into
`unclassified`. `DesignParameterDerivative` contains `dA/dp`, frame/reference
orbit derivatives, and any entrance/exit boundary terms. An edge-angle report
is not accepted if it contains only the local quadrupole profile and omits this
chain.

The report can rank intervals for one matrix entry such as `R12`, `R34`, or a
coupling entry such as `R13`, or for a weighted matrix objective. This makes the
cause of linear coupling observable by path position, normal/skew multipole,
longitudinal/fringe field, element, and field source. `NonlinearMapReport`
extends this to `T_ijk`, `U_ijkl`, direct local nonlinearities, and ordered
two-region cascades without hiding the raw coefficients behind one quality
score.

## pybind11 Surface

The public Python package is `radia.beam`. C++ bindings may remain internal in
`_radia_pybind`, but Python classes must be thin owners of the native objects.

The implemented HCurl vector-potential entry is:

```python
result = radia.beam.propagate_grid_function_linear_map(
    solved_a,
    lengths_m,
    reference_positions_m,
    reference_tangents,
    magnetic_rigidity_t_m,
    sample_radius_m=1e-3,
    initial_horizontal=[1.0, 0.0, 0.0],
    names=region_names,
)

k1 = result["normal_gradient_per_m2"]
k1s = result["skew_gradient_per_m2"]
R = result["R"]
```

The caller supplies one exact reference position and tangent per constant
segment. Reference-orbit construction is intentionally not hidden inside this
first adapter.

The implemented inspectable tracking foundation is:

```python
from radia import beam

species = beam.ParticleSpecies.proton()
reference = beam.ReferenceParticle.from_kinetic_energy_ev(species, 220e6)
field = beam.UniformField(magnetic_t=[0.0, 0.0, 0.7])
initial_state = beam.CartesianState(
    position_m=[0.0, 0.0, 0.0],
    kinetic_momentum_kg_m_s=[reference.momentum_kg_m_s, 0.0, 0.0],
)
sample = field.evaluate([0.0, 0.0, 1.2])
equation = beam.LorentzEquation(species, field, independent="time")
rhs = equation.rhs(0.0, initial_state)

stepper = beam.Boris2()
one_step = stepper.step(equation, 0.0, initial_state, 1e-12)
tracker = beam.Tracker(equation, stepper)
plan = beam.TrackPlan()
plan.start = 0.0
plan.stop = 1e-8
plan.maximum_step = 1e-11
trajectory = tracker.track(initial_state, plan)
```

The following remains the target composition once general field/lattice and
closed-orbit objects are implemented:

```python
field = beam.SumField([
    beam.GridFunctionField(solved_a),  # HCurl A; native curl(A) is B
    beam.TransformedField(error_b, alignment),
])

problem = beam.ClosedOrbitProblem(tracker, lattice, frame)
residual = problem.residual(initial_guess)
closed = problem.solve(initial_guess, return_history=True)

analysis = beam.TransferAnalyzer(field, reference).analyze(trajectory)
quad = analysis.normal_quadrupole_profile
skew = analysis.skew_quadrupole_profile
r13 = analysis.local_strength_sensitivity(
    "skew_quadrupole", matrix_entry=(0, 2))
edge = analysis.design_sensitivity(
    edge_angle_parameter, objective=target_transfer_matrix)

nonlinear = analysis.analyze_nonlinear(
    trajectory, maximum_order=3,
    amplitude_box=[2e-3, 5e-3], delta_range=[-0.01, 0.01])
edge_nonlinear = nonlinear.design_sensitivity(
    edge_angle_parameter,
    penalize=["T", "U", "chromatic", "resonance"])
```

NumPy arrays are accepted and returned for batches, maps, and trajectories.
Object properties expose native metadata and immutable configuration. Python
callbacks are permitted only in an explicitly named diagnostic adapter; they
are not a production or per-step backend.

## Standalone MEX Surface

MATLAB uses the same `radia_mex` command namespace. Expensive immutable fields
and lattices use checked `uint64` handles. State, steps, residuals, stencils,
maps, and diagnostics are ordinary MATLAB structs and numeric arrays.

Required commands:

```text
beam.field.from_grid_function
beam.field.info
beam.field.sample
beam.field.destroy

beam.lattice.create
beam.lattice.info
beam.lattice.destroy

beam.equation.rhs
beam.step
beam.track

beam.closed_orbit.residual
beam.closed_orbit.jacobian
beam.closed_orbit.solve

beam.transfer_map.stencil
beam.transfer_map.track_stencil
beam.transfer_map.fit
beam.transfer_map.build

beam.linearization.local_generator
beam.multipoles.profile
beam.transfer.propagate_variational
beam.transfer.from_grid_function
beam.transfer.local_strength_sensitivity
beam.transfer.design_sensitivity
beam.transfer.interval_summary
beam.nonlinear_map.build
beam.nonlinear_map.validate
beam.nonlinear_attribution.analyze
beam.nonlinear_attribution.region_pairs
beam.normal_form.analyze

beam.optics.analyze
beam.attribution.analyze
```

The MEX registry validates handle type, generation, liveness, and dependency
ownership. `beam.field.sample`, `beam.equation.rhs`, and `beam.step` are
independently callable and tested. A MATLAB class such as `radia.beam.Field`
owns a handle and destroys it in `delete`; it is convenience, not a second
implementation.

The standalone MEX ABI is the debugging and parity surface. A future Simulink
S-Function calls these kernels but does not become their only entry point.

## No EarlyTimes Format Compatibility

The migration intentionally does not parse or emit `et.inp`,
`latticeParameter.inp`, `table.dat`, `input-TOSCA.dat`, Fortran namelists,
fixed-column lattice records, numbered Fortran units, or EarlyTimes `.xy`
files. It also does not expose the old six-value cylindrical state array.

Applications construct `ParticleSpecies`, `ReferenceParticle`, `Field`,
`Element`, `Lattice`, and solver option objects directly. Python uses typed
classes and NumPy arrays; MATLAB uses typed structs, numeric arrays, and native
handle classes. General-purpose importers may be added for maintained public
formats, but they are independent of the EarlyTimes migration.

## Existing GridFunction And `.sol` Boundary

This section records the field handoff used by the analyzer. It is supporting
infrastructure, not the EarlyTimes migration's principal beam-dynamics result.
The transfer-attribution API consumes a live field object and is independent of
whether that object was solved in-process or restored from a checked bundle.

The preferred path is an in-process native vector-potential `GridFunction`:

```python
p = 5
A = ngsolve.GridFunction(ngsolve.HCurl(mesh, order=p), name="A")
# Project only the gauge-constrained HDiv-MMM A CoefficientFunction.
maps = project_earlytimes_grid_function_maps(A_coefficient, mesh)
A = maps.vector_potential
# Keep B_coefficient as the independent HDiv-MMM direct-B RK source.
```

The `GridFunction` retains its `FESpace`, and the space retains its mesh.
Python therefore does not pass a raw DoF vector, separately declared order, or
sampled regular field to the tracker.  The independent B-map tracker evaluates
the original B CoefficientFunction and never derives it from A.  An optional
HDiv projection is requested explicitly with
`project_magnetic_flux_density=True`; it is a projection check, not the B
source of truth.

For a process or language boundary, `.sol` is supported as part of a checked
solution bundle, not as a standalone field file:

```text
beam_field/
  mesh.vol
  A.sol
  space.json
  bundle.json
```

`space.json` records the NGSolve space type, order, dimension, vector dimension,
real/complex type, `definedon`, relevant flags, and expected DoF count.
`bundle.json` records the NGSolve and Radia versions, vector-potential unit
(`T m`), derived magnetic-field unit (`T`), mesh and file digests, coordinate
convention, field name, and bundle schema version.

Python uses native NGSolve operations:

```python
mesh = ngsolve.Mesh("mesh.vol")
space = make_space_from_checked_spec(mesh, "space.json")  # HCurl(order=p)
A = ngsolve.GridFunction(space, name="A")
A.Load("A.sol")
field = radia.beam.GridFunctionField(A)
```

MATLAB uses the same native C++ path without Python:

```matlab
space = radia.ngsolve.FESpace.create("mesh.vol", "hcurl", order);
A = radia.ngsolve.GridFunction.fromFESpace(space, Name="A");
A.load("A.sol", "space.json", "bundle.json");
result = radia.beam.propagateGridFunctionMultipoleMap( ...
    A, lengthsM, positionsM, tangents, magneticRigidityTM);
```

The required MEX additions are:

```text
ngsolve.grid_function.save
ngsolve.grid_function.load
ngsolve.grid_function.evaluate
beam.field.from_grid_function
```

`save` and `load` call NGSolve's native `GridFunction::Save` and
`GridFunction::Load`; MATLAB never parses `.sol`. Loading fails before assigning
data when the mesh digest, FESpace signature, DoF count, value type, vector
dimension, units, or supported NGSolve version differs.

The beam commands use the mesh already retained by the native `GridFunction`
handle and evaluate its HCurl derivative directly. They do not require a
second mesh path and do not copy the DoF vector into MATLAB.

NGSolve `.sol` data does not by itself define the mesh or finite-element space.
This is why a bare `.sol` path is not accepted by `GridFunctionField`. For
multilevel or refined meshes, the bundle producer must freeze the final mesh
and exact space before saving; arbitrary reconstruction from an unrelated
`.vol` is rejected. A native NGSolve object archive may later be added for
cases where the refinement hierarchy itself must be preserved.

## Definition Of Complete Migration

"Complete" means all canonical, reachable numerical behavior linked by the
EarlyTimes `src/Makefile` is classified in a checked coverage manifest and
either:

1. implemented by a named C++ component;
2. delegated to a named NGSolve or Radia component with the same numerical
   meaning; or
3. classified as unreachable historical/debug code with evidence; or
4. classified as retired format-specific parsing, printing, or file control.

The required feature families are:

| EarlyTimes capability family | C++ owner |
|---|---|
| initialization and configuration state | typed constructors and option objects |
| ideal, hard-edge, analytic, and solved 3-D fields | analytic fields and `NGSolveGridFunctionField` |
| field scaling, placement, and errors | field expression tree |
| global/local coordinate transforms | explicit frame conversions |
| time, azimuth, and segment tracking | `Equation`, `Stepper`, `Tracker` |
| RK4 and second-order Lie tracking | named `Stepper` implementations |
| reference and closed orbit | `OneTurnMap`, `ClosedOrbitProblem` |
| multipole expansion and fitting | `TransferMapBuilder`, multipoles |
| distributed linear/nonlinear map attribution | `TransferAnalyzer`, `NonlinearMapReport` |
| beta, dispersion, tune, FFT, resonance | optics functions |
| COD and correction | correction problem and solver |
| injection, extraction, septum, kicker | elements, thin operators, events |
| generated distributions and scattering | explicit source models with RNG |
| monitor values and plots | result value objects and maintained exporters |

Files named `*_old`, `*_obsolete`, `.g`, prebuilt objects/binaries, disabled
`if (.false.)` experiments, and working-directory outputs are inventoried but
are not copied into production source. A still-reachable call wins over its
filename and must be migrated.

## Verification Gates

1. **Core invariants:** zero-field drift, uniform magnetic circle, momentum
   conservation, charge-sign reversal, frame round-trip, and event location.
2. **Distributed linear optics:** drift, normal quadrupole, skew quadrupole,
   sector bend, solenoid, and analytic edge/fringe fixtures reproduce their
   closed-form local generators and endpoint maps.
3. **Attribution identity:** `R(sf,s) R(s,s0) = R(sf,s0)` at every station;
   integrated sensitivity densities match direct Frechet derivatives and small
   design perturbations. Region summaries sum to the reported derivative.
4. **Frame objectivity:** parallel-transport and explicitly rotated frames
   transform normal/skew coefficients and maps consistently without creating
   false physical coupling.
5. **Nonlinear analytic maps:** ideal sextupole and octupole fixtures reproduce
   the expected `T` and `U` coefficients; translated reference orbits reproduce
   analytic feed-down; longitudinal and momentum terms retain the declared
   factorial and index conventions.
6. **Nonlinear attribution:** station densities reconstruct `T`; direct and
   ordered cascade terms reconstruct `U`; order-by-order Poisson-bracket tests
   verify the canonical map through the requested order.
7. **Validity domain:** direct tracking over the declared amplitude/momentum
   box agrees with the truncated map within its reported tolerance and fails
   closed outside the accepted domain.
8. **Field parity:** analytic magnetic fields and native NGSolve `curl(A)` from
   `HCurl(order=p)` GridFunctions, including FESpace-class/order metadata and
   checked `.sol` bundle round trips in Python and MATLAB. Explicit direct-B
   mode remains an independent comparison path.
9. **Integrator parity:** recorded canonical EarlyTimes trajectories for RK4
   and Lie routes, plus convergence-order tests.
10. **Orbit and map parity:** closed orbit, first/second-order map coefficients,
   tune, beta, dispersion, and footprint from preserved EarlyTimes outputs.
11. **Transparency:** each high-level result must be reproducible from its
   returned field samples, steps, stencil rays, or iteration history.
12. **Binding parity:** C++, Python, and MEX produce the same values, statuses,
   errors, metadata, and storage conventions.
13. **Lifecycle:** stale, wrong-type, destroyed, and cross-kernel MEX handles fail
   loudly; repeated create/track/destroy cycles leak no native state.
14. **Performance:** warmed C++, pybind11, and MEX timings are measured on the
   same field and trajectory after numerical parity passes.

Heavy historical cases belong in `validation_test/early_times/`. Fast analytic
and API tests belong in `tests/`. Preserved reference outputs are converted to
small, unattributed golden fixtures rather than copying an entire old working
directory.

## Migration Order

1. Freeze this API and generate the EarlyTimes routine coverage manifest.
2. Implement value types, fields, frames, Lorentz equations, RK4, Boris, and
   trajectory records in C++.
3. Implement local multipole fitting, canonical derivatives `A/B/C`,
   forward/reverse variational maps `R/T/U`, station/interval sensitivity
   densities, and ordered nonlinear cascade attribution. Cross-check first
   against the existing planar combined-function Frechet chain and analytic
   sextupole/octupole fixtures.
4. Add pybind11 and standalone MEX parity for `field.sample`, local
   linearization, variational propagation, and transfer sensitivity before
   exposing a high-level tracking orchestration command.
5. Port the ideal FFAG and hard-edge algorithms and connect solved NGSolve
   `GridFunction` fields.
6. Port closed orbit, reference-orbit storage, transfer maps, and linear optics.
7. Port nonlinear optics, tune/FFT/resonance, COD, injection, and extraction.
8. Add actual CAD/topology parameter derivatives, objective-ranked edge/fringe
   recommendations, and Simulink composition.
9. Mark complete only when the coverage manifest has no unclassified reachable
   routine and the C++/Python/MEX gates pass on one commit.

## Source Collaboration And Attribution Gate

Y. Ishi confirmed that using the calculation and tracking design presents no
licensing problem. AI-assisted conversion or mechanical translation of the
EarlyTimes implementation requires consultation with the original author and
should use the latest source supplied by the author.

Until that consultation and source handoff are complete, the archived source is
used only to inventory capabilities and design the public API. Do not translate
the archived routines into production C++, and do not copy prebuilt binaries,
object files, historical snapshots, parsers, or working outputs into Radia.

Before implementation conversion begins, record:

1. the latest-source delivery date, version or commit, and checksum;
2. the agreed AI-assisted conversion and public-release scope;
3. the attribution text requested by the original author; and
4. the numerical behavior selected for preservation or deliberate retirement.

Public papers that use the EarlyTimes-derived tracking design cite Sheehy et
al., IPAC'15 paper MOPJE077, including its publisher and DOI. The citation is an
attribution requirement for publications; it does not make the obsolete
EarlyTimes file formats part of the Radia API.

## References

- [NGSolve CoefficientFunctions](https://docu.ngsolve.org/ngs24/tutorials/02_coefficientfunction.html)
- [NGSolve custom CoefficientFunction add-on pattern](https://docu.ngsolve.org/v6.2.2202-18-gcd039cbb4/mylittlengs/1_Basic/1_myFEM/coefficient.html)
- EarlyTimes ET User Guide, NFG-KM-0030, 2003-10-14 (lab archive)
- S. L. Sheehy et al., "Progress on Simulation of Fixed Field Alternating
  Gradient Accelerators," in *Proc. IPAC'15*, Richmond, VA, USA, 2015,
  pp. 495-498, paper MOPJE077,
  [doi:10.18429/JACoW-IPAC2015-MOPJE077](https://doi.org/10.18429/JACoW-IPAC2015-MOPJE077).
