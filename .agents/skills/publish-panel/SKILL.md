---
name: publish-panel
description: Verify and publish the independently released Radia Simulink application library through its package-owned gates.
---

# Publish application blocks

Run `simulink-app-health` on the tracked library and MATLAB-only samples under
`packages/radia-simulink`. Preserve the five production application contracts,
native MEX lifecycle and required numerical/parity checks. Use the supported
MATLAB tooling for model edits, save/reopen and diagnostics.

Build and verify the package-owned distribution before publication. Simulink is
independent of Radia's four-case release wrapper; Radia wrapper completion does
not certify MATLAB/MEX/SLX artifacts. Cubit's toolbar remains separately owned.
Never add Qt dependencies to normal Radia Python.
