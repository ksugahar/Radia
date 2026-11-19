# MCP capability packs

Capability packs combine existing domain tools in one stdio process. The domain
modules remain the implementation owners. Tool names, schemas, annotations,
structured results, images, prompts and static resources retain their original
contracts. Calls use the SDK's public APIs and retain the domain call-log and
reload wrappers. There are no child MCP transports or generic tool-execution
proxies in the client-visible tool list.

## Served packs

Every pack runs as `radia-mcp serve <pack>`.

| Pack | Profiles (in addition to `all`) |
| --- | --- |
| `radia-mcp-analysis` | `field`, `verification`, `theory` |
| `radia-mcp-optimize` | `search`, `topology`, `coil` |
| `radia-mcp-learning` | - (GNN, PINN, data assimilation) |
| `radia-mcp-ih` | - (IH, force, MOR) |
| `radia-mcp-maglev` | - (MagLev, electromagnet, force) |
| `radia-mcp-electromagnet` | - (electromagnet, accelerator, materials) |
| `radia-mcp-motor` | - (motor, electromagnet, materials) |
| `radia-mcp-pcb` | - (PCB, PEEC, litz/transmission) |
| `radia-mcp-acoustic` | `production`, `education` |
| `radia-mcp-postprocess` | `blender` (Gmsh stays in `gmsh-mcp`) |
| `radia-mcp-rna` | - (RNA/MEC, magnetic materials and hysteresis) |

The default is `all`. Use a narrow profile when one exists, e.g.
`radia-mcp serve radia-mcp-analysis --profile field`. Selection happens before domain imports. Discovery is fixed at startup;
restart after adding/removing tools or changing a profile. Existing body-only
reload behavior belongs to each domain. A duplicate tool, prompt or resource
identifier aborts startup rather than choosing an implementation silently.

`capability_pack_status` reports the active profile, source, version and exact
tool owners. Source-domain status tools still describe the domain, not the whole
pack. The acoustic production and educational routes remain explicitly named;
neither substitutes for the other. Application packs group existing workflows
and verification tools; they do not claim a newly implemented autonomous solver.

## Client migration

After installing the wheel, replace old client entries with
`radia-mcp serve <pack>` and optional `--profile` arguments, keyed by the pack
name. For example, replace the separate optimization server entries with one
`radia-mcp-optimize` entry. Reconnect the client
after editing its configuration. Client configuration changes are separate from
package installation and are not performed by this implementation.

Cubit, build123d, Gmsh, official MATLAB/Optuna, Radia MATLAB and Mathematica
remain independent process/installation boundaries. Choose only the capability
packs needed for a task; combining servers alone does not reduce schema size.

## Verification

Run `pytest tests/test_capability_packs.py` from the package, plus the affected
catalog tests. Tests compare complete live schemas and cover real stdio
initialization, discovery, domain calls, prompts, resource reads and errors.
Every new command supports `--selftest` (profile import and discovery only).
Wheel checks must run outside the checkout without `PYTHONPATH`; use a separate
environment so another session's editable installation is not redirected.

New pack members must have no server-owned startup/shutdown lifecycle. Resource
templates currently fail explicitly until a routing adapter is supplied. Do not
add session-owning external servers to the manifest. Keep numerical kernels and
independent validation methods in their owning modules.
