# radia-mcp examples

Runnable smoke / demo scripts that exercise the radia-mcp MCP tools
end-to-end. Useful as both regression tests and as worked-example
references when authoring AI prompts.

| Script | Demonstrates |
|---|---|
| [`cadquery_to_cubit_hex_demo.py`](cadquery_to_cubit_hex_demo.py) | CadQuery L-bracket → STEP → `cubit_mesh_auto` ladder → live Cubit GUI. End-to-end one-call pipeline. |
| [`cubit_exec_safely_demo.py`](cubit_exec_safely_demo.py) | `cubit_exec_safely` workflow: build user state → AI clean recipe (commits) → AI broken recipe (blocked at batch, GUI untouched) → checkpoint listing. |
| [`gmsh_post_v22_to_v41_demo.py`](gmsh_post_v22_to_v41_demo.py) | `gmsh_post_*`: inspect v2.2 Cubit export → validate (fails on version) → convert to v4.1 → validate (passes) → write `$NodeData` view → re-inspect. |

## Running

Each script is self-contained and assumes radia-mcp is installed
(`pip install radia-mcp[full]`) and Coreform Cubit is on the
default install path. Run with:

```bash
python examples/cadquery_to_cubit_hex_demo.py
```

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
