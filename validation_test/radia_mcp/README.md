# radia-mcp Numerical Validation

This directory contains solver-backed numerical validation for the
`radia_mcp.radia_ngsolve` APIs. These checks generate Netgen/NGSolve meshes,
run FEM or BEM solves, compare formulations, or exercise convergence and
frequency sweeps. They are intentionally outside the fast package test suite.
Gmsh figure/raster checks that generate STEP geometry, probe volume fields, or
render numerical images live here for the same reason; package tests retain
only solver-free option and public-signature contracts.

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
validation budget or requires hibino, with mdx as an idle-CI fallback.

## MCP server surface performance

`benchmark_mcp_server_surface.py` compares an untouched baseline checkout with
the candidate `radia-mcp` source. It records cold-import latency, direct tool
count, serialized `tools/list` schema size, and `full`-profile compatibility in
`results_mcp_server_surface_20260902.json`. This is explicit performance
validation and is intentionally excluded from fast CI. The comparison preserves
every supported baseline tool in the full profile while requiring the explicitly
listed retired standalone-panel tools to remain absent; any other missing tool is
a validation failure.
