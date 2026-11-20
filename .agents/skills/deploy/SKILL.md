---
name: deploy
description: Route Radia deployment to its owning release workflow and verify the selected live runtime.
---

# Deploy

For Radia, use `release-quad`: a release waits for the `cpu`,
`linux-python-gpu`, `windows-python-gpu` and `matlab` lane receipts of one SHA,
and the solver runtime is maintained on the validation role only. The CPU-only
build host prepares CPU candidates and receives no GPU runtime. Follow
[editable operation and isolated QA](../../../docs/operations/editable-deployment.md).
Host-specific skill names and runtime mappings are private configuration.

Use `verify-deploy` to confirm fresh imports, intended source and native artifact
identity. Keep normal editable operation separate from isolated exact-wheel
acceptance. Resolve explicit host/interpreter/source bindings before mutation;
preserve existing environments, other packages and active jobs.

MCP follows its [independent runtime policy](../../../packages/radia-mcp/docs/operations/mcp-runtime-policy.md).
Cubit follows `release-cubit-mesh-export` on its licensed build role only.
Simulink follows `simulink-app-health` and its package-owned release gates.
A Radia release does not install or certify these independent packages.

Use the reviewed owning controller. Do not bypass failed target checks, rebuild
native artifacts on a consumer, or silently substitute a CPU/GPU backend.
