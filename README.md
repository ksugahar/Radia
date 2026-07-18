# Radia - AI-Native Electromagnetic CAE

[![CI](https://github.com/ksugahar/Radia/actions/workflows/build-test.yml/badge.svg)](https://github.com/ksugahar/Radia/actions/workflows/build-test.yml)
[![Policy Lint](https://github.com/ksugahar/Radia/actions/workflows/policy-lint.yml/badge.svg)](https://github.com/ksugahar/Radia/actions/workflows/policy-lint.yml)
[![Python 3.12](https://img.shields.io/badge/python-3.12-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-see%20LICENSE-blue.svg)](LICENSE)

**AI designs. Radia provides the engineering platform.**

Radia is a programmable electromagnetic CAE platform built on
[NGSolve](https://ngsolve.org/). It connects Python, MCP, Jupyter, CAD,
meshing, electromagnetic analysis, optimization, and visualization into one
engineering workflow.

Radia is not another all-in-one solver and it is not a replacement for
NGSolve. NGSolve remains the numerical foundation. Radia contributes the
electromagnetic methods, open-boundary models, application workflows, and
AI-facing interfaces that are needed to turn a numerical backend into a
practical design platform.

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
            +--> Electromagnetic design and analysis
            |       Radia methods + NGSolve / ngsolve.bem
            |
            +--> Optimization, validation, and visualization
                    Python, Jupyter, Gmsh, Netgen, result artifacts

An LLM can write and execute the Python workflow through MCP. A human can
inspect the same workflow in a notebook or application panel. The interfaces
are different; the engineering model and validation artifacts are shared.

## Architecture

Radia is organized as three layers.

| Layer | Responsibility |
| :--- | :--- |
| **Application** | Magnet design, Hodograph, VIM, Eddy, Stream Function, induction heating, MagLev, motors, and other concrete workflows |
| **Platform** | Python APIs, MCP servers, Jupyter workbenches, build123d, Cubit, Gmsh, Netgen, validation, and result artifacts |
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
6. Jupyter, Gmsh, or Netgen presents the durable result to a human engineer.

Humans remain in the loop for assumptions, physical interpretation, and
release decisions. AI automation is valuable because the workflow is
executable and inspectable, not because it removes engineering judgment.

### radia-mcp

The [radia-mcp package](packages/radia-mcp/) provides MCP servers and
knowledge tools for the Radia ecosystem. It covers Radia and NGSolve
workflows as well as Cubit, Gmsh, build123d, PEEC, optimization, analytical
references, panel workbenches, and validation.

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
runtime.

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

## Jupyter workbenches

The canonical human-facing panels are notebook workbenches. They wrap the
same headless calculation scripts used by automation and produce durable
run.log and result.json artifacts.

    python -m jupyter lab src/radia/panels/notebooks/radia_ih.ipynb
    python -m jupyter lab src/radia/panels/notebooks/radia_em.ipynb
    python -m jupyter lab src/radia/panels/notebooks/radia_pcb.ipynb
    python -m jupyter lab src/radia/panels/notebooks/radia_streamfunction.ipynb

The Cubit toolbar is a separate, Cubit-embedded integration surface. Normal
Radia Python workflows do not install or depend on Cubit's private PySide
runtime.

## Documentation

- [Documentation index](docs/README.md)
- [HDiv-VIM](docs/hdiv_vim/README.md)
- [Clebsch-Hodograph](docs/clebsch_hodograph/README.md)
- [Eddy-current methods](docs/solver/EDDY_CURRENT_METHODS.md)
- [Stream Function](docs/stream_function/README.md)
- [PEEC integration](docs/peec_integration/README.md)
- [Induction heating](docs/induction_heating/README.md)
- [Electric machines](docs/electric_machine/README.md)
- [API reference](docs/api/API_REFERENCE.md)
- [Panel development](docs/panels/ADDING_NEW_PANEL.md)
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
