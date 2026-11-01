# Curved Kelvin artifacts: MCP discovery and evidence

Named periodic export is useful when two curved boundary patches must represent
the same translated geometry, rather than merely share matching corner vertices.
The exporter owns that geometry; NGSolve owns the finite-element spaces and maps.
Radia MCP explains the consumer checks through its existing Kelvin lookup tools.
It does not maintain a second mesher or a patched Netgen runtime.

This page is an introduction and discovery map. Current operating instructions,
supported representations and failure semantics come from the selected live
server's status, catalog and knowledge responses. Do not copy a generated tool
inventory or a past status response into a second procedural manual.

## Use the existing package-owned servers

| Purpose | Existing entry point | Read-only discovery |
|---|---|---|
| Cubit CAD, neutral capture and mesh export | `mcp-server-cubit` from `cubit-mesh-export` | `tools/list`, `cubit_status`, `cubit_topics` |
| Kelvin formulation and solver-artifact guidance | `mcp-server-radia-ngsolve` from `radia-mcp` | `tools/list`, the registered status tool, `kelvin_transformation` and `kelvin_identify_post_hoc` |

These are independent distributions. Cubit-only operation does not require
`radia-mcp`. No additional server, wrapper package or version bump is needed to
make an inventory. Use the actual selected interpreter and client configuration;
an entry point found on a different PATH does not identify a running client.
Follow the [maintenance procedure](maintenance.md) for client setup and source
identity checks. Reconnect only the owned client if its loaded source is stale;
preserve other clients, CAD sessions and jobs.

## Find the producer and consumer contracts

Start with the exporter's [MCP guide](../../cubit-mesh-export/src/cubit_mesh_export/mcp/README.md)
and [named periodic export contract](../../cubit-mesh-export/README.md#explicit-named-periodic-export).
The supported named-curved route uses one translational `kt` PERIODIC group and
polynomial q1/q2/q3 tetrahedra with linear connectivity. It does not qualify
arbitrary orders, rotations, independent high-order-node meshes or chained
multi-direction seam ownership.

Through the existing Radia server, ask for:

```text
kelvin_transformation(topic="curved_geometry")
kelvin_transformation(topic="identify")
kelvin_identify_post_hoc(topic="curved_geometry")
```

These calls select documentation only. The existing post-hoc tool also has a
`verify` branch that loads a mesh and adds in-memory identifications; its
tool-wide annotation therefore does not promise read-only behavior. Do not use
that execution branch as a discovery, high-order geometry or persistence check.

The lookup covers master-derived slave coefficients, q3 orientation parity,
shared-entity ownership, complete pair/patch coverage and independent map checks.
The final `.vol` requires its verified `.vol.json` companion record. Pair count,
curve-order integers or a successful build alone are insufficient evidence.

For imported Cubit high-order archives, full-Netgen Load/Save is not an automatic
post-hoc repair. Persistence findings depend on the exact route and runtime;
native OCC reconstruction is a different route. Keep element transformations
alive during mapped-point use. A translation-consistent curved seam need not
lie on the exact analytical inversion sphere. Geometry validity, HCurl trace
continuity, field accuracy, p-convergence and release qualification are distinct.
The live knowledge response owns the detailed checks and applicability limits.

## Keep neutral capture separate from external adapters

The public [model-capture contract](../../cubit-mesh-export/docs/model_capture.md)
records source CAD and explicit units/identities without meshing. It is not
arbitrary exact CAD translation, and materials or physical boundary conditions
are not guessed from names. A separately owned external adapter must qualify its
supported reconstruction and saved/reopened output. Exporter publication does
not certify that adapter or make its private runtime part of this public server.
Private engineering models and comparison evidence do not belong in this guide.

## Close the user workflow, not only publication

A complete handoff distinguishes these records:

1. The package's reviewed source, artifact identity and publication checks.
2. The actual client's selected interpreter, imported source/version and native
   artifact identity where applicable, obtained through its live status route.
3. A harmless `tools/list` and knowledge lookup on that client, with the expected
   producer/consumer limitations visible. Static source tests do not establish
   that a deployed client has loaded the change.
4. Linked public documentation, affected checks and independent review on the
   integrated source. A mesh or solver execution needs its own artifact ownership
   and validation evidence; a read-only lookup does not certify numerics.

Keep operational receipts with the task's private evidence. Do not call a
published tag, an installed wheel, a private draft or another client's response
proof that all four are complete. Preserve models and evidence during cleanup;
remove a worktree only after its unique changes and active consumers are resolved
through the assigned management workflow.
