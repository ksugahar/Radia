---
name: acoustic-fembem
description: Work on the separate readable MATLAB acoustic FEM-BEM education solver. Use for .vol meshes, P1 FEM/BEM, convolution quadrature, or MATLAB-vs-NGSolve teaching validation; production acoustics remains in radia.acoustics.
---

# MATLAB Acoustic FEM-BEM

## Stack Contract

Use the official MATLAB MCP Server from MathWorks as the MATLAB execution
substrate:

- `detect_matlab_toolboxes` for MATLAB/toolbox availability.
- `check_matlab_code` before accepting nontrivial `.m` code.
- `evaluate_matlab_code` for short probes.
- `run_matlab_file` for scripts/examples.
- `run_matlab_test_file` for MATLAB unit tests.

Do not replace this layer with a Sugahara-specific MATLAB bridge. The separate
education solver adds domain knowledge above the official server: `.vol` IO, P1 FEM/BEM
assembly checks, convolution quadrature conventions, Gypsilab-like API shape,
and NGSolve / `ngsolve.bem` cross-validation.

This package does not implement the production `radia-acoustic` MATLAB API.
Production scattering, FSI, and CQ live in `radia.acoustics`; their common
MATLAB boundary belongs to `radia_mcp.matlab`.

Radia owns the custom-tool extension contract. Locate and inspect it with:

- `acoustic_fembem_extension_contract` from `mcp-server-acoustic-fembem`;
- `radia-acoustic-fembem-mcp-config --contract-only` from a shell.

Generate official-server arguments with
`acoustic_fembem_server_config` or `radia-acoustic-fembem-mcp-config`. The
companion GPL MATLAB teaching solver remains outside the BSD `radia-mcp`
wheel; set `ACOUSTIC_FEMBEM_ROOT` or pass `--project-root` when those MATLAB
entry points are needed.

If the official MATLAB MCP tools are not available in the current agent session,
say so and use the nearest safe fallback (`matlab -batch` or foreground SSH to
mdx for long jobs). Do not silently pretend a MATLAB run happened.

## Default Workflow

1. Identify the MATLAB project root and expected entry point:
   - library code: `+package/`, `@class/`, or ordinary function folders;
   - tests: `tests/` with `matlab.unittest`;
   - examples or apps: separate from reusable solver kernels.

2. Probe the environment:
   - MATLAB version and toolboxes via `detect_matlab_toolboxes`;
   - required paths with a short `evaluate_matlab_code` probe;
   - `.vol` / mesh availability before solver work.

3. Implement MATLAB in small, testable units:
   - pure assembly and quadrature functions first;
   - thin scripts/apps second;
   - no path mutation except in explicit setup helpers.

4. Validate:
   - `check_matlab_code` on changed `.m` files;
   - `run_matlab_test_file` for unit tests;
   - one independent numerical cross-check whenever feasible.

5. Report exact numerical settings:
   mesh, element order, quadrature, frequency/time grid, CQ weights, solver
   tolerance, conditioning/rank observations, and validation tolerance.

## Domain Gates

- `.vol` is the durable mesh interchange format. Do not invent a parallel mesh
  text format for MATLAB unless the task explicitly requires it.
- P1 FEM/BEM means geometry, DOF ordering, mass/stiffness signs, boundary normal
  conventions, and singular/near-singular BEM quadrature must be checked before
  trusting a plot.
- Gypsilab-style APIs should be familiar but not copied blindly. Preserve the
  useful shape (mesh, domain/boundary, finite element space, integral operator)
  while naming lab-specific acoustic or Radia concepts explicitly.
- Convolution quadrature work must pin the Laplace-domain convention, time step,
  generating function, and inverse transform normalization in tests.
- Cross-validation priority:
  1. analytical small geometry;
  2. NGSolve / `ngsolve.bem` Python path;
  3. independent MATLAB implementation;
  4. regression golden from a previous verified run.

## Promotion Rule

Keep experimental MATLAB notebooks/scripts in `C:\temp` until the method is
stable. Promote reusable code to `src/`, durable validation to
`validation_test/<topic>/`, human-facing explanation to `docs/<topic>/`, and
agent-facing rules to this skill or a `radia_mcp.<topic>` server.

When a MATLAB workflow becomes a repeated agent operation, add a thin
`radia_mcp` extension tool that prepares inputs, calls the official MATLAB MCP
execution layer or a documented MATLAB entry point, and returns structured
validation results. The extension should not become a second MATLAB runtime
manager.
