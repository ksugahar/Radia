# Start with an engineering question

**Discover a capability, then operate CAE with AI.**
Radia is electromagnetic CAE you can operate through an AI assistant using MCP.
You decide the analysis purpose, physical conditions and acceptance criteria;
the assistant can invoke supported CAD, mesh, calculation and result operations.
Python and MATLAB/Simulink interfaces remain available, with NGSolve as the
finite-element foundation.

[Application examples](APPLICATION_GUIDE.md) · [Accuracy and validation](VALIDATION_GUIDE.md) ·
[Install Radia](../README.md#quick-start) · [All documentation](README.md)

## What AI operation changes

- **Less tool-specific operation to learn.** Describe the requested analysis
  operation instead of memorizing every command or menu sequence.
- **Reuse existing inputs.** Start from supported CAD, mesh or model data.
  Check transferred units, materials, sources and boundaries before solving;
  importing a shape alone does not transfer an entire physical model.
- **Connect the workflow.** Use MCP-accessible CAD and mesh operations, then
  execute the selected analysis and inspect its artifacts. Mesh controls remain
  explicit rather than being hidden behind a generic automatic mesh.

The current capability is AI-driven CAE operation. The user determines the
analysis strategy and judges the result; Radia does not claim autonomous design.

## Choose your first result

The [six selected applications](APPLICATION_GUIDE.md#browse-results) show coils,
motor ripple, local heating losses, magnetic lift and particle trajectories.
For detailed numerical work, continue to the
[validation study index](../validation_test/STUDY_INDEX.md).

| Your goal | First stop | What you can inspect |
| :--- | :--- | :--- |
| Calculate a magnet's field | [First field calculation](../README.md#quick-start) | A magnetic-field value in tesla, without an air mesh |
| Learn the field-model workflow | [NGSolve integration notebook](ngsolve_integration/integration_basics.ipynb) | Saved geometry, mesh and field results |
| Design a coil for a target field | [Coil-design case](APPLICATION_GUIDE.md#design-a-coil-for-a-target-field) | Current distribution, field target and winding workflow |
| Understand conductor loss and circuit models | [Conductor case](APPLICATION_GUIDE.md#connect-conductors-to-circuits) | Skin-effect examples and circuit-model workflow |
| Study electromagnetic heating | [Heating case](APPLICATION_GUIDE.md#follow-loss-into-temperature) | Source transfer and transient temperature verification |
| Build a dynamic motor model | [Motor case](APPLICATION_GUIDE.md#inspect-a-reduced-motor-model) | Angle-dependent flux and torque comparisons |
| Judge whether a result is trustworthy | [Validation guide](VALIDATION_GUIDE.md) | Reference solutions, measured quantities and scope |

## Choose how you want to work

| You want to… | Begin with | First useful result |
| :--- | :--- | :--- |
| Collaborate on CAD and mesh in Cubit | [Cubit mesh workflow](cubit_mesh_export/README.md) | Geometry, mesh controls, labels and checked solver inputs |
| Study small executable samples | [tests/](../tests/README.md) | API calls, fixtures, assertions and expected behavior |
| Study physics and application samples | [validation_test/](../validation_test/README.md) | Model setup, numerical checks and scoped evidence |
| Inspect results before installing | [Visual application gallery](APPLICATION_GUIDE.md#browse-results) | A geometry or response plot, its conditions and the saved notebook |
| Write a small Python calculation | [Analytical quick start](../README.md#quick-start) | A magnetic-field value in tesla, without CAD or mesh setup |
| Use NGSolve objects from MATLAB | [Mesh, space and matrix walkthrough](../matlab/README.md#persistent-mesh-space-form-and-matrix-handles) | Selected native objects and operations available from MATLAB |
| Connect an application to Simulink | [Radia library guide](../matlab/README.md) | The application block's inputs, outputs and execution contract |
| Ask an AI assistant to operate CAE tools | [Python and MCP setup](../README.md#python-and-mcp) | Your analysis task carried out through supported tools, with inspectable artifacts |

Follow the linked setup for the chosen interface. A MATLAB MEX operation, a
triggered batch application and a native Simulink time-step block have different
execution requirements; see the [backend map](../matlab/python_api_parity_manifest.json).
These links do not imply that every workflow is available as a standalone
MATLAB Runtime application.

## Explore capabilities

**Read the model conditions, inspect the field, understand the supported conclusion.**

1. Use the [small magnet example](../README.md#quick-start) to establish input
   units and read one computed field value.
2. Open the [NGSolve integration notebook](ngsolve_integration/integration_basics.ipynb)
   to relate model geometry, mesh and field. Inspect the saved views before rerunning.
3. Choose a [worked application](APPLICATION_GUIDE.md#read-a-result-before-choosing-a-method)
   and read its conditions, measured output and interpretation.
4. Apply the [validation guide](VALIDATION_GUIDE.md) to your changed geometry or
   excitation before drawing a design conclusion.

## From feature introduction to operation

1. **Inspect a saved result.** Open an application notebook before installing
   dependencies. Saved outputs show what the example computed;
   interactive views may require a compatible notebook frontend.
2. **Choose an operating route.** Use [MCP](../README.md#python-and-mcp), or
   collaborate with AI using [Cubit](cubit_mesh_export/README.md) and
   [Simulink](../matlab/README.md). Notebook execution serves the demonstration;
   it is not a required step in operating your model. AI Cubit execution is
   headless; geometry, journals, meshes and checks support the handoff to human work.
3. **Check the inputs.** Record geometry, SI units, material law, excitation,
   boundary conditions and the quantity you need. Mesh workflows require the
   labels and domains specified by the selected application.
4. **Check the result.** Compare an appropriate reference, refine the mesh or
   time step, and retain the result artifacts. The [validation guide](VALIDATION_GUIDE.md)
   explains the difference between solver convergence and physical accuracy.

## Before you install

The [installation guide](../README.md#installation) owns supported platforms,
dependencies and package choices. Radia, radia-mcp and the other distributions
have separate release histories; installing one does not install every interface.
The [MATLAB backend map](../matlab/python_api_parity_manifest.json) identifies
native operations and explicit Python batch fallbacks.

## Common questions

**Where are the executable samples?**
[tests/](../tests/README.md) contains small usage and contract examples;
[validation_test/](../validation_test/README.md) contains physics and application
cases. The repository uses the singular name `validation_test/`. Read each case
with its fixtures, dependencies and assertions before adapting it.

**Do I operate Radia by editing documentation notebooks?**
Notebooks introduce features and explain saved calculations. Normal operation
uses MCP or human–AI collaboration through Cubit and Simulink. You may rerun a
demonstration to study it within its stated scope and dependencies.

**Does AI decide what to design or which analysis is appropriate?**
You define the purpose, assumptions and analysis plan, and assess the engineering
conclusions. AI can help carry out supported operations through MCP. Optimization
is an available calculation workflow when you specify its objectives and constraints;
it is not a promise of autonomous design.

**Do I need MATLAB or a CAD license for the first example?**
The analytical Python field example does not require either. MATLAB/Simulink
and licensed CAD workflows have their own requirements. Mesh-based workflows
also have an explicit build123d/Netgen route; follow the chosen application's inputs.

**Does Radia replace NGSolve?**
NGSolve supplies the finite-element foundation. Radia adds electromagnetic
operators, open-boundary methods, application coupling and engineering workflows.

**Is every example production-qualified?**
A notebook demonstrates its saved configuration. Numerical validation and
release acceptance are separate evidence. Each application page should identify
its limitations; the [project status](../README.md#project-status) describes the overall maturity.

**Where do I go when a calculation fails?**
Open a [bug report](https://github.com/ksugahar/Radia/issues) with the smallest
shareable model, package versions, command or block configuration, expected
quantity, and relevant result/error record. Remove credentials and private model
data. See [contribution guidance](../CONTRIBUTING.md) and the private
[security reporting process](../SECURITY.md).

**How do I keep up with changes?**
Read the [changelog](../CHANGELOG.md) and [release notes](https://github.com/ksugahar/Radia/releases).
For method-level details, continue to the [documentation index](README.md).
