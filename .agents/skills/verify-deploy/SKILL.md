---
name: verify-deploy
description: Verify Radia normal editable runtime identity separately from isolated exact-wheel acceptance.
---

# Verify deployed Radia

Use the reviewed `tools/release_quad.py verify-editable` with explicit private
host, absolute interpreter, source and exact-commit bindings. Both normal
Windows roles use editable Radia. Follow
[editable operation and isolated QA](../../../docs/operations/editable-deployment.md).

In a fresh selected interpreter, verify editable `direct_url.json`, actual
Radia import origin/version, exact source commit, native payload hashes and ABI,
metadata-pinned dependencies, `pip check` and the recorded Radia source intent.
A mismatch fails; do not repoint implicitly or manufacture an intent record.

The QA interpreter must be a separate isolated venv holding the exact accepted
wheel. Verify wheel/installed RECORD hashes and execute the required adopted
CPU/GPU acceptance there. Normal editable identity is not wheel qualification.
Compilation/import alone is not numerical or architecture qualification.

MCP, Cubit's embedded toolbar and MATLAB/Simulink have independent runtime and
verification contracts. Verify each requested consumer separately with its
owning skill. Preserve other packages, environments and active jobs.
