# radia-mcp examples

Runnable smoke / demo scripts that exercise the radia-mcp MCP tools
end-to-end. Useful as both regression tests and as worked-example
references when authoring AI prompts.

| Script | Demonstrates |
|---|---|
| [`cadquery_to_cubit_hex_demo.py`](cadquery_to_cubit_hex_demo.py) | CadQuery L-bracket → STEP → `cubit_mesh_auto` ladder → live Cubit GUI. End-to-end one-call pipeline. |
| [`cubit_exec_safely_demo.py`](cubit_exec_safely_demo.py) | `cubit_exec_safely` workflow: build user state → AI clean recipe (commits) → AI broken recipe (blocked at batch, GUI untouched) → checkpoint listing. |
| [`gmsh_post_v22_to_v41_demo.py`](gmsh_post_v22_to_v41_demo.py) | `gmsh_post_*`: inspect v2.2 Cubit export → validate (fails on version) → convert to v4.1 → validate (passes) → write `$NodeData` view → re-inspect. |
| [`dtn_spectrum_coarse_mesh_demo.py`](dtn_spectrum_coarse_mesh_demo.py) | Why open boundaries are accurate on coarse meshes, MEASURED two ways (NGSolve, no Cubit). **Part A** `exterior_dtn_spectrum`: the BEM matrix Λ_h = V⁻¹(−½M+K) on a sphere at two mesh sizes → eigenvalues vs the ladder −(n+1)/R → the three spectral facts (low modes accurate coarse, error ordered by degree, refinement widens the band). **Part B** `kelvin_dtn_eigenvalue`: the Kelvin closure's effective DtN by volume FEM — mode n inverts to a degree-n polynomial, exact iff FEM order≥n (the dominant dipole → linear → order-1 coarse accurate). **Part C** `kelvin_openbc_error_vs_exterior_mesh`: accuracy vs the EXTERIOR (Kelvin ball) mesh size with the interior FIXED — ISOLATES the Kelvin open-boundary error from the interior FEM error (swap only the Γ operator). The open-BC error converges ~×4/level as the exterior refines but stays 45–709× below the fixed ~5% interior FEM error, so a coarse exterior mesh already suffices. Kameari's exterior-refinement accuracy check as a quantified, separated error budget — see MCP tool `dtn_coarse_mesh`. |

## Running

Each script is self-contained and assumes radia-mcp is installed
(`pip install radia-mcp[full]`) and Coreform Cubit is on the
default install path. Run with:

```bash
python examples/cadquery_to_cubit_hex_demo.py
```

## Heavy Validation

Long FEM/physics cross-validations are not part of the default
radia-mcp pytest gate. They are marked `xval` and should be run
explicitly when validating a research workflow:

```bash
python -m pytest tests/ -m xval
```

When a validation becomes primarily explanatory or publication-facing,
move it from `tests/` into `examples/` as a runnable script with a short
README entry. Keep `tests/` focused on fast MCP/API contracts and small
numerical invariants.

The Cubit-touching demos kill any existing `coreform_cubit.exe`
process at start (so they leave a clean slate); the live Cubit GUI
window is left running for inspection at the end.

## Conventions

- All test artifacts are written into a fresh `tempfile.mkdtemp(...)`
  directory, **not** into the user's working directory or the repo.
- Scripts are idempotent: re-running them should produce the same
  pass/fail outcome.
- Output uses ASCII / UTF-8 so the demos work in cp932 console
  environments without crashing.
