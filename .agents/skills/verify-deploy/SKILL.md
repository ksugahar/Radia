---
name: verify-deploy
description: Verify that the Radia validation runtime holds the exact verified release wheel.
---

# Verify deployed Radia

The validation role runs verified wheels only; editable installs exist only on
the build host. Follow
[validation-role verified-wheel deployment](../../../docs/operations/validation-wheel-deployment.md).

Use the reviewed `tools/release_quad.py phase9` with explicit private host,
hostname and absolute interpreter bindings. In a fresh selected interpreter it
refuses an editable `direct_url.json` and checks the Radia import origin,
release version, ABI pins and native payload digest. `phase8` additionally
binds every installed RECORD hash to the hash-verified candidate wheel.
A mismatch fails; repair forward by deploying the verified release wheel with
`phase8`, never an editable source.

Execute the required adopted CPU/GPU acceptance on that wheel runtime.
Compilation/import alone is not numerical or architecture qualification.

MCP, Cubit's embedded toolbar and MATLAB/Simulink have independent runtime and
verification contracts. Verify each requested consumer separately with its
owning skill. Preserve other packages, environments and active jobs.
