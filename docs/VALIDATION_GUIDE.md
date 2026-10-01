# How to judge a Radia result

**A useful result identifies what was compared, under which conditions, and with what error.**
Use this guide to move from a feature introduction to its numerical evidence.
Executable samples live in [tests/](../tests/README.md) for small usage and
contract checks and [validation_test/](../validation_test/README.md) for physics
and application cases. Read their setup and assertions together. Documentation
notebooks explain selected results; they do not replace those sample collections
or the MCP and Cubit/Simulink operating routes.

[Start here](START_HERE.md) · [Application examples](APPLICATION_GUIDE.md) · [All documentation](README.md)

## Choose evidence for the quantity you need

| Engineering question | Reference or check | Evidence entry point |
| :--- | :--- | :--- |
| Does current diffuse correctly into a conductor? | Closed-form round-wire resistance, cylinder penetration, half-space skin effect and lamination loss | [Diffusion suite](../validation_test/eddy_current_analytical_validation/README.md#diffusion-closed-form-suite) and [recorded results](../validation_test/eddy_current_analytical_validation/diffusion_closed_form_results.json) |
| Is a force or torque calculation consistent? | Lorentz force, coenergy/virtual work, Maxwell-stress balance and symmetry | [Force and torque corpus](../validation_test/force_validation/README.md) |
| Is the thermal discretization and source transfer correct? | Manufactured transient solution and conservative transfer | [Thermal verification notebook](induction_heating/axisymmetric_p2_thermal.ipynb) and [scope](induction_heating/README.md) |
| Does a motor reduced model generalize beyond its sampled angles? | Interlaced holdout angles and independent torque evaluations | [Motor-ROM study](electric_machine/angle_periodic_motor_rom.ipynb) and [qualification notes](electric_machine/README.md) |
| Do Python and MATLAB perform the same supported operation? | Explicit interface comparison cases | [Parity record](../validation_test/ngsolve_matlab_parity/README.md) and [backend classifications](../matlab/python_api_parity_manifest.json) |

These links point to existing evidence, not a claim that every check was rerun
for this guide. Read the artifact's date, source/runtime identity, inputs,
reference, measured error and threshold before applying a result elsewhere.
Historical failures and restricted applicability remain part of the evidence.

## Three different questions

**Did the equations converge?** A small algebraic residual says the discrete
system was solved to its specified tolerance. It does not measure model error.

**Did the discretization resolve the quantity?** Mesh, polynomial-order,
quadrature and time-step studies reveal sensitivity. A smooth field picture
alone does not establish convergence, especially for force and loss.

**Does the model represent the intended physics?** Compare an analytic solution,
a manufactured solution, an independently derived formulation, or an appropriate
public reference. Record material laws, excitation, boundary conditions and the
measurement region. Agreement between two implementations sharing the same
assumption cannot test that assumption.

## Reading a numerical claim

- **Quantity and units:** field vector, integrated loss, force, torque, temperature or another stated observable.
- **Error definition:** relative or absolute; RMS, maximum or integral; sample set or physical region.
- **Operating range:** geometry, frequency, current, temperature and material conditions covered.
- **Reference:** exact solution, independent discretization, measured data or a stored regression record.
- **Reproduction:** generator/driver, inputs, versions and saved numerical result.
- **Limitations:** unresolved cases, extrapolation, unsupported geometry and missing acceptance steps.

For reduced models, report full-model discretization error separately from
reduction error. Validate inputs or parameter points excluded from training.
For speed claims, match accuracy and report offline preparation as well as
online cost under comparable execution conditions.

## Demonstration, numerical evidence and release acceptance

| Artifact | What it establishes | What it does not establish by itself |
| :--- | :--- | :--- |
| Executed notebook in `docs/` | A visible calculation for saved inputs | Universal accuracy or production readiness |
| Numerical record in `validation_test/` | A measured check with its declared scope | Coverage outside the tested conditions |
| Interface/backend map | How a capability is implemented and exposed | Numerical equivalence for every use |
| Release acceptance | The specified package/runtime passed its required gates | Validity of a new engineering model |

For a new application, start with a nearby analytical case, reproduce its result,
then change one physical assumption at a time. If a result is surprising, retain
the input and error record and use the [support guidance](START_HERE.md#common-questions).
