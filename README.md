# Radia - AI-Native Electromagnetic CAE

[![CI](https://github.com/ksugahar/Radia/actions/workflows/build-test.yml/badge.svg)](https://github.com/ksugahar/Radia/actions/workflows/build-test.yml)
[![Policy Lint](https://github.com/ksugahar/Radia/actions/workflows/policy-lint.yml/badge.svg)](https://github.com/ksugahar/Radia/actions/workflows/policy-lint.yml)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-see%20LICENSE-blue.svg)](LICENSE)

**AI designs. Radia provides the engineering platform.**

Radia is a programmable electromagnetic CAE platform built on
[NGSolve](https://ngsolve.org/). It connects Python, MCP, Simulink, CAD,
meshing, electromagnetic analysis, optimization, result-bearing notebooks,
and visualization into one engineering workflow.

Radia is not another all-in-one solver and it is not a replacement for
NGSolve. NGSolve remains the numerical foundation. Radia contributes the
electromagnetic methods, open-boundary models, application workflows, and
AI-facing interfaces that are needed to turn a numerical backend into a
practical design platform.

> **Use Radia through its MCP servers.** This README is the entry point, not a
> duplicate operating manual. Ask the relevant `radia-mcp` server for the
> current workflow, supported contract, machine-specific prerequisites, and
> validation steps before running a domain operation. The MCP response is
> designed to be executable by an agent and routes MATLAB execution through
> MathWorks' official MATLAB MCP Server. Start with the
> [radia-mcp tool catalog](packages/radia-mcp/README.md).

## Highlights

- **Couple motion, control, circuits, and electromagnetic force in Simulink.**
  Simulink owns the time-domain motion and controller model, while Radia
  supplies electromagnetic force/field models and LTspice supplies the drive
  circuit. Vector signals and interval state handoff close the loop between
  position, velocity, current, voltage, and force.
- **Use existing LTspice circuits as the electrical subsystem.** Run real
  device models as a sampled-data circuit plant without Simscape Electrical,
  with hidden execution and capacitor-voltage/inductor-current handoff between
  intervals. The hysteretic plant block iterates LTspice current and Radia
  magnetic flux within each interval, returning hysteretic back EMF to the
  power-electronic circuit before committing either state.
- **Optimize the coupled machine, circuit, and controller together.** Run
  parallel Simulink/LTspice trials with MATLAB-native Optuna-style studies,
  multiple objectives, and live Pareto monitoring. Native NGSolve/Cubit sheet
  topology can use the same trial lifecycle through `SheetMetalRunner`.
- **Turn a SPICE netlist into an editable LTspice schematic.** Convert `.cir`
  or `.net` files to automatically laid-out `.asc` schematics from Python or
  MATLAB, then verify connectivity by converting back with LTspice.
- **Use the same engineering APIs from Python and MATLAB.** Radia keeps its C++
  numerical kernels and Python circuit parser as single sources of truth while
  exposing checked MATLAB and Simulink adapters.
- **Use NGSolve from MATLAB through a native MEX bridge.** Persistent checked
  handles expose selected NGSolve meshes, finite-element spaces, coefficient
  and grid functions, forms, vectors, sparse matrices, and HCurl workflows to
  MATLAB without launching Python for each operation. The bridge preserves
  NGSolve's ownership of finite-element mathematics and is cross-checked
  against the Python path.

```text
             position / velocity
        +------------------------------+
        |                              v
Simulink motion <--- force --- Radia / NGSolve electromagnetic model
        ^                              ^
        |                              | current
        +--- controller ---> LTspice drive circuit
                    gate / voltage command
```

```matlab
% Netlist -> editable LTspice schematic, with round-trip connectivity check
radia.ltspice.netlistToSchematic("converter.cir", ...
    OutputFile="converter.asc", ValidateRoundTrip=true);

% Use the LTspice drive circuit inside a Simulink motion/control model
radia.simulink.buildLTspiceBlock("current_controller", ...
    Netlist="converter.cir", ...
    InputNames=["gate"; "reference"], ...
    OutputTraces=["I(L1)"; "V(out)"], ...
    SampleTime_s=1e-4);
```

Register the packaged blocks once after adding Radia's `matlab` directory to
the MATLAB path:

```matlab
radia.simulink.buildLibrary
sl_refresh_customizations
```

The Simulink Library Browser then shows one **Radia** library containing the
Electromagnet, PCB PEEC, Motor, Stream Function, Induction Heating, LTspice,
Optuna, and Sheet Metal Optimization blocks. For the machine-readable setup and compatibility contract, ask the
`mcp-server-radia-matlab` tool `matlab_simulink_library_contract`.

## The central idea

> NGSolve owns the numerical foundation.
> Radia adds the missing electromagnetic engineering methods.
> AI and humans use the same programmable workflow.

The platform is designed for a complete engineering loop:

    Natural-language intent
            |
            v
    Python / MCP workflow
            |
            +--> CAD and mesh generation
            |       build123d, Cubit, Netgen, Gmsh
            |
            +--> Electromechanical co-simulation
            |       Simulink motion/control + LTspice circuits
            |       + Radia/NGSolve electromagnetic force
            |
            +--> Electromagnetic design and analysis
            |       Radia methods + NGSolve / ngsolve.bem
            |
            +--> Optimization, validation, and visualization
                    Python, Jupyter, Gmsh, Netgen, result artifacts

An LLM can write and execute the Python workflow through MCP. A human operates
the application through Simulink and inspects durable results in docs notebooks
or native visualization tools. Notebook workbenches are retired for every
application, including IH.
The interfaces differ; the engineering model and validation artifacts are shared.

## Architecture

Radia is organized as three layers.

| Layer | Responsibility |
| :--- | :--- |
| **Application** | Magnet design, Hodograph, VIM, Eddy, Stream Function, induction heating, MagLev, motors, and other concrete workflows |
| **Platform** | Python APIs, MCP servers, Simulink application blocks, docs notebooks, build123d, Cubit, Gmsh, Netgen, validation, and result artifacts |
| **Numerical** | NGSolve finite elements and ngsolve.bem for spaces, transformations, quadrature, weak forms, BEM operators, and linear algebra |

The boundary between these layers matters. Radia should extend NGSolve at the
physics and application layers, while continuing to use NGSolve's public
abstractions for finite-element plumbing.

## What Radia owns

Radia focuses on electromagnetic capabilities that are not provided by a
general-purpose finite-element backend alone:

- analytical magnetic and source fields for open regions;
- magnetic-material methods based on volume integral and charge-Gram ideas;
- surface impedance, PEEC, and low-frequency eddy-current workflows;
- topology-aware coil and conductor design;
- model reduction for repeated electromagnetic solves;
- orchestration that lets AI agents use CAD, meshing, analysis, and validation
  as one workflow.

## What Radia delegates

Radia deliberately integrates strong existing tools instead of rebuilding
them:

- **NGSolve**: finite-element spaces, Piola maps, curved geometry,
  orientations, quadrature, weak-form assembly, GridFunctions, and the
  numerical solve;
- **ngsolve.bem**: boundary-element formulations and surface operators;
- **Netgen**: mesh generation and notebook visualization;
- **Coreform Cubit**: CAD-driven mesh generation and the Cubit export
  interface;
- **build123d**: programmable CAD construction;
- **Gmsh**: durable mesh and field-result visualization;
- **NumPy, SciPy, MKL, and proven linear-algebra libraries**: numerical
  building blocks.

This is not only a packaging preference. It keeps Radia's development effort
focused on electromagnetic methods rather than duplicating mature numerical
infrastructure.

## Core technologies

### Electromechanical co-simulation with LTspice and Simulink

Radia uses Simulink as the dynamic-system orchestrator: mechanical position and
velocity update the electromagnetic model, electromagnetic force updates the
motion model, and the LTspice drive circuit exchanges current, voltage, and
control signals with both. This motion-circuit-field loop does not require
Simscape or Simscape Electrical. LTspice support is part of that coupled
workflow rather than the top-level purpose by itself:

- bidirectional LTspice conversion between SPICE netlists (`.cir`, `.net`)
  and editable LTspice schematics (`.asc`);
- a Python-native converter in `radia-spice-lab`, exposed to MATLAB through a
  thin wrapper so parsing, symbol dictionaries, and automatic layout have one
  implementation;
- hidden LTspice execution on Windows, including SSH-driven runs on the LAB
  and 100号機 development hosts;
- real and complex RAW import, stepped analyses, `.noise`, and transient FFT
  APIs;
- a distributable **Radia / LTspice Circuit** block in the Simulink Library
  Browser with vector inputs and outputs;
- sampled-data co-simulation that advances LTspice one interval at a time and
  hands saved node voltages and inductor currents to the next interval;
- a transactional **Hysteretic LTspice Plant** block that waveform-iterates
  circuit current, vector Play/Energy hysteresis, flux linkage, and back EMF;
  failed iterations roll back both circuit and material state, while converged
  steps output current, flux density, flux linkage, back EMF, Maxwell force,
  and hysteresis energy;
- MATLAB-native Optuna-style single- and multi-objective studies, parallel
  Simulink/LTspice trials, native sheet-metal topology trials, and live Pareto
  monitoring with standard Simulink Scope and XY Graph blocks;
- analytic-adjoint continuous optimization with bounded MMA or gradient-only
  SQP, including a direct native HCurl Eddy-Bubble material-topology adapter.

```matlab
% Convert a netlist to an editable schematic and verify connectivity by
% converting it back with LTspice.
schematic = radia.ltspice.netlistToSchematic("plant.cir", ...
    OutputFile="plant.asc", ValidateRoundTrip=true);

% Add a state-handoff LTspice plant with two control inputs and two outputs.
radia.simulink.buildLTspiceBlock("controller_model", ...
    Netlist="plant.cir", ...
    InputNames=["gate"; "reference"], ...
    OutputTraces=["V(out)"; "I(L1)"], ...
    SampleTime_s=1e-4);
```

The Python circuit converter remains the source of truth; MATLAB and Simulink
delegate to it and to the same LTspice execution/RAW contracts. See the
[MATLAB integration guide](matlab/README.md) and
[radia-spice-lab documentation](packages/radia-spice-lab/README.md).

### Hodograph

Hodograph methods turn selected nonlinear magnetic design problems into a
linearized design problem in a transformed coordinate space. They provide a
direct route to magnetic flux-line and pole-face design, especially for
accelerator and precision magnet workflows.

- [Clebsch-Hodograph documentation](docs/clebsch_hodograph/README.md)

### VIM

The VIM is Radia's successor to magnetic-moment methods. It is not an ELF
compatibility layer. Its purpose is to evolve open-boundary magnetic analysis
around superposition, linearity, NGSolve-compatible spaces, and practical
low-frequency applications.

The current HDiv-VIM path uses NGSolve meshes and finite-element spaces while
Radia supplies the electromagnetic charge-Gram and open-boundary operators.
This makes the method suitable for soft magnetic materials, nonlinear
magnetization workflows, and coupled magnetic applications.

- [HDiv-VIM documentation](docs/hdiv_vim/README.md)

### Eddy

Radia's Eddy framework is a high-order edge-element framework for
eddy-current problems. The Eddyable concept identifies basis functions and
reduced models that are specialized for the electromagnetic response of a
particular problem, reducing the cost of repeated solves without hiding the
underlying NGSolve formulation.

- [Eddy-current methods](docs/solver/EDDY_CURRENT_METHODS.md)

### Stream Function

The Stream Function layer supports topology-aware coil and conductor design.
It connects field objectives, regularization, contour extraction, and
manufacturable single-stroke coil paths to the broader optimization workflow.

- [Stream Function documentation](docs/stream_function/README.md)

## Application layer

The application layer is where the methods become engineering tools.

| Workflow | Typical use |
| :--- | :--- |
| **Radia Magnet** | Permanent magnets, coils, accelerator magnets, and open-space field design |
| **Radia VIM** | Magnetic materials, soft iron, nonlinear demagnetization, and coupled magnetics |
| **Radia Eddy** | Low-frequency eddy currents, SIBC, ESIM, shielding, and reduced transient response |
| **Radia Stream Function** | Coil topology, winding design, field shaping, and optimization |
| **PEEC and circuit workflows** | Inductance, resistance, coupling, skin/proximity effects, and SPICE-compatible extraction |
| **Induction heating and motors** | Application-specific workflows built from the common Radia and NGSolve layers |

Representative application domains include magnetic levitation, wireless power
transfer, induction heating, accelerator magnets, motors, printed-circuit
conductors, and other open-space magnetic systems.

## Respect for NGSolve

NGSolve is not a dependency that Radia wraps casually. It is the numerical
language and foundation that Radia builds upon.

Radia follows these principles:

1. Use NGSolve spaces, forms, GridFunctions, and mapped evaluation APIs for
   finite-element work.
2. Let NGSolve own element orientation, local-to-global transformations,
   Piola mappings, curved geometry, quadrature, and weak-form assembly.
3. Add Radia-specific physics around those abstractions rather than
   reimplementing finite-element plumbing in Python.
4. Keep independent analytic, integral, or reduced routes where they improve
   validation and physical insight.
5. Prefer a clear NGSolve workflow over a Radia-specific parallel vocabulary
   when NGSolve already provides the right abstraction.

In short: Radia extends NGSolve; it does not compete with it.

## Physical scope

Radia targets magneto-quasi-static and Darwin-regime electromagnetic
problems. Radia's interaction kernels are Laplace kernels, with surface
impedance and skin depth handling frequency-dependent conductor physics.
Radia is not a full-wave Helmholtz solver and does not aim to replace
full-wave tools for radiation-dominated problems.

The combination is useful when:

- the air region is large or effectively unbounded;
- a magnet or coil moves without wanting to remesh the surrounding air;
- conductor skin depth would make a volume mesh impractical;
- the magnetic source is best represented analytically;
- many parameter variations or optimization steps are required.

## AI-native workflow

AI is a first-class user of the platform.

The intended workflow is:

1. An LLM turns an engineering request into a parameterized Python model.
2. MCP tools construct CAD, generate or inspect meshes, and select the
   appropriate Radia and NGSolve workflow.
3. A headless calculation script runs the solve and writes run.log and
   result.json.
4. Analytic references, mesh checks, and independent formulations validate
   the result.
5. Optimization varies the design while preserving the model and its
   provenance.
6. Simulink provides the operating surface; docs notebooks, Gmsh, or Netgen
   present the durable result to a human engineer.

Humans remain in the loop for assumptions, physical interpretation, and
release decisions. AI automation is valuable because the workflow is
executable and inspectable, not because it removes engineering judgment.

### radia-mcp

The [radia-mcp package](packages/radia-mcp/) provides MCP servers and
knowledge tools for the Radia ecosystem. It covers Radia and NGSolve
workflows as well as Cubit, Gmsh, build123d, PEEC, optimization, analytical
references, Simulink application contracts, and validation.

Install it separately when an MCP client is available:

    python -m pip install radia-mcp

Start with the Radia metadata/catalog server to discover the available
domain servers and tools. The package README contains client-specific
configuration examples.

## Integration routes

| Need | Preferred route |
| :--- | :--- |
| CAD construction | build123d or Coreform Cubit |
| Mesh generation | Netgen/NGSolve or Cubit export to Netgen .vol |
| FEM | NGSolve |
| BEM and surface operators | ngsolve.bem |
| Open-boundary source fields | Radia analytical field APIs |
| Magnetic materials | HDiv-VIM coupled to NGSolve |
| Eddy currents | NGSolve HCurl workflows, ngsolve.bem, SIBC, ESIM, and reduced models |
| Field and mesh visualization | Netgen WebGUI and Gmsh/GmshPostExport |
| Circuit extraction | Radia PEEC and optional radia-spice-lab |
| Optimization | Python scientific stack, parameter sweeps, and optimization libraries |

For Cubit-to-NGSolve workflows, the .vol file is the process boundary:
Cubit produces the mesh and the NGSolve/Radia process consumes it. This keeps
the Cubit Python 3.10 runtime separate from the Radia/NGSolve Python 3.12
runtime. Run `check-vol` after export and before solver or Simulink
initialization; it checks NGSolve structure, curved mappings, labels, and an
optional Cubit CAD sidecar. Production modes add a versioned strict label
contract. Physical material constants remain explicit DesignSpec/configuration
values and are never inferred from mesh labels.

## Quick start

### Supported development target

| Component | Target |
| :--- | :--- |
| Operating system | Windows 10/11 or Windows Server |
| Python | 3.12 |
| NGSolve / Netgen | 6.2.2604 |
| Coreform Cubit | 2025.12, optional |
| Native build | MSVC + Intel MKL |

### Install the Python package

For the core Python API and NGSolve integration:

    python -m pip install radia

For the optional Cubit mesh-export plugin:

    python -m pip install "radia[cubit]"
    cubit-plugin-install --verify-only

For AI-assisted workflows:

    python -m pip install radia-mcp

Pin all package versions together for a reproducible lab deployment. The
release workflow validates the compatibility of radia,
cubit-mesh-export, and radia-mcp across the supported machines.

### A minimal analytic field

Radia uses SI units. Magnetization is specified in A/m and the field is
returned in tesla.

    import numpy as np
    import radia as rad

    mu0 = 4.0 * np.pi * 1e-7
    remanence = 1.2  # tesla

    magnet = rad.ObjRecMag(
        [0.0, 0.0, 0.0],
        [0.01, 0.01, 0.01],
        [0.0, 0.0, remanence / mu0],
    )

    field = rad.Fld(magnet, "b", [0.0, 0.0, 0.02])
    print("B [T] =", field)
    rad.UtiDelAll()

The source field is evaluated analytically. There is no air mesh in this
minimal example.

### NGSolve source coupling

Radia fields can be supplied to NGSolve as native coefficient functions:

    import radia as rad

    source_B = rad.RadiaField(magnet, "b")
    source_A = rad.RadiaField(magnet, "a")

    # Use source_B or source_A in the appropriate NGSolve
    # GridFunction, BilinearForm, or LinearForm workflow.

The NGSolve mesh, finite-element space, mapped evaluation, and assembly stay
in NGSolve. Radia supplies the electromagnetic source term.

### Native NGSolve MEX bridge for MATLAB

Radia also makes selected NGSolve workflows directly usable from MATLAB. The
native MEX gateway keeps checked `uint64` handles for `Mesh`, `FESpace`,
`CoefficientFunction`, `GridFunction`, `BilinearForm`, `LinearForm`, `Vector`,
and `Matrix` objects. MATLAB can assemble finite-element operators, inspect
sparse matrices in NGSolve's global degree-of-freedom ordering, and run
HCurl-based reduced and topology workflows without starting a Python process
for each operation.

This bridge is valuable in its own right: it lets MATLAB and Simulink users
build on NGSolve's mature finite-element foundation while sharing Radia's C++
electromagnetic kernels and artifact contracts. It is not a fork or a
reimplementation of NGSolve, and it does not claim to duplicate the complete
NGSolve Python API. NGSolve still owns spaces, orientation, transformations,
Piola maps, curved geometry, quadrature, assembly, and field evaluation.

The MEX contracts are tested against the corresponding Python/NGSolve results,
including complex state/adjoint conventions, matrix layout, and HCurl
multifrequency topology gradients. Acoustic soft, rigid, fluid, and elastic
sphere references plus convolution-quadrature primitives also share one C++
implementation through pybind11 and MEX. Full acoustic CQ-BEM and FSI retain
NGSolve's Python object model and use an explicit Python-DLL fallback from
MATLAB rather than duplicating finite-element plumbing. Axisymmetric Henrotte
FEM has begun the same promotion: Q1 stiffness and conductivity-mass element
matrices share one C++ implementation across the NGSolve BFI, pybind11, and
MEX, while NGSolve continues to own spaces and global assembly. See the
[MATLAB integration guide](matlab/README.md) for the supported surface, build
instructions, and runtime requirements.

When a field objective exposes an analytic sensitivity, MATLAB can pass the
native MEX adjoint directly to `radia.topopt.optimizeAdjoint` and select MMA or
SQP. `radia.topopt.optimizeHCurlActivationAdjoint` already closes the native
multifrequency HCurl Joule adjoint with a material-volume constraint. TPE and
CMA-ES remain the outer/global choices for discrete and conditional designs;
MMA/SQP perform the continuous inner refinement without numerical-gradient
fallbacks.

## Simulink application blocks

The canonical human-facing interfaces are masked blocks in the single Radia
Simulink library: Electromagnet, PCB PEEC, Motor, Stream Function, and
Induction Heating. Electromagnet, PCB PEEC, Motor, and Stream Function delegate
to the same `DesignSpec` and headless calculation used by Python/MCP. Induction
Heating is the native exception: its Eddy and Thermal blocks are C/C++ MEX
S-Functions with no Python fallback in the simulation loop.

The four batch blocks write `run.log` and `result.json`; field-producing runs
also write a checked GMSH `.msh v4.1` artifact in the run directory. The IH
preview currently exposes its field vectors directly in Simulink and does not
claim that batch artifact contract.

    addpath("matlab")
    radia.setup()
    radia.simulink.buildLibrary()

The four batch application blocks launch the validated Python CLI only on a
rising trigger; they do not start Python every simulation step. Promoting those
blocks to MEX/ROM remains optional and requires parity and long-run testing.
The native NGSolve MEX bridge above is a separately supported platform
capability, not merely a future block backend.

In MATLAB, `radia.simulink.openIH()` opens `radia_ih.slx`. The first native IH
release is explicitly a preview runtime for checked, preassembled Eddy and
Thermal operators. It proves the S-Function lifecycle, current and rotation
inputs, conservative temperature transport, and closed thermal feedback. The
remaining Cubit `.vol` to PEEC/BEM-A/BIM/FEM operator-assembly boundary is not
claimed complete by that preview. LUT-only and lumped state-space models are
not shipped as Radia IH interfaces.

The Cubit toolbar is a separate, Cubit-embedded integration surface. Normal
Radia Python workflows do not install or depend on Cubit's private PySide
runtime.

## Documentation

- [Documentation index](docs/README.md)
- [Executable Radia <-> NGSolve example with saved WebGUI scenes](docs/ngsolve_integration/integration_basics.ipynb)
- [HDiv-VIM](docs/hdiv_vim/README.md)
- [Clebsch-Hodograph](docs/clebsch_hodograph/README.md)
- [Eddy-current methods](docs/solver/EDDY_CURRENT_METHODS.md)
- [Stream Function](docs/stream_function/README.md)
- [PEEC integration](docs/peec_integration/README.md)
- [Induction heating](docs/induction_heating/README.md)
- [Electric machines](docs/electric_machine/README.md)
- [API reference](docs/api/API_REFERENCE.md)
- [Simulink application-block development](docs/panels/ADDING_NEW_PANEL.md)
- [Cubit mesh export](docs/cubit_mesh_export/README.md)
- [Build from source](BUILD.md)
- [MCP package and tool catalog](packages/radia-mcp/README.md)

## Development

Clone the repository and install the editable packages:

    git clone https://github.com/ksugahar/Radia.git
    cd Radia
    python -m pip install -e .
    python -m pip install -e packages/cubit-mesh-export
    python -m pip install -e packages/radia-mcp

The native extension is built with MSVC, Intel MKL, CMake, and Ninja:

    pwsh -NoProfile -ExecutionPolicy Bypass -File .\Build.ps1 -RadiaOnly

Run focused tests while developing, then broaden the test and validation
scope according to the change:

    python -m pytest -q tests/test_vim_eddy_hybrid.py
    python -m pytest -q

Fast regression tests live under tests/. Important numerical checks and
research-grade sweeps live under validation_test/. Result-bearing method
notebooks belong under docs/.

## Project status

Radia is an active research and engineering platform. The technical
foundation is being strengthened before broadening the application surface.

Current priorities:

1. strengthen the platform and its executable validation;
2. complete and document the NGSolve-compatible VIM workflow;
3. establish robust Hodograph design workflows;
4. improve the Eddy and model-reduction framework;
5. consolidate user documentation as the architecture stabilizes.

The repository should lead with analytic solutions, physical methods, and
reproducible workflows. Local validation provenance and machine-specific
deployment notes stay in the development environment rather than public
user-facing material.

## Heritage and license

Radia originates from the magnetostatics work developed at the European
Synchrotron Radiation Facility. The current project extends that heritage
with NGSolve integration, open-boundary engineering methods, high-order
workflows, MCP interfaces, and AI-oriented automation.

Radia contains components with different license terms. See [LICENSE](LICENSE)
for the complete terms, including the Radia core, HACApK, and sparsesolv
components.
