# radia-mcp Numerical Validation

This directory contains solver-backed numerical validation for the
`radia_mcp.radia_ngsolve` APIs. These checks generate Netgen/NGSolve meshes,
run FEM or BEM solves, compare formulations, or exercise convergence and
frequency sweeps. They are intentionally outside the fast package test suite.

Run the fast API and MCP contract tests with:

```powershell
python -m pytest packages/radia-mcp/tests -m "not xval and not slow"
```

Run this numerical validation lane from the repository root with:

```powershell
python -m pytest validation_test/radia_mcp -m "not slow and not compute_host"
python -m pytest validation_test/radia_mcp
```

The root `validation_test/conftest.py` automatically applies the `validation`
marker. Add `slow` or `compute_host` only when a case exceeds the normal manual
validation budget or requires an idle mdx/hibino compute host.
