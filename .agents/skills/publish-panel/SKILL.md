---
name: publish-panel
description: Publish Radia's Simulink application library after verifying its five application blocks and tracked MATLAB-only samples. Use simulink-app-health before deploy/release-qud; Cubit's toolbar remains a separate deployment surface.
---

# Publish Application Blocks

1. Run `simulink-app-health`.
2. Regenerate and inspect `matlab/radia_simulink_library.slx`.
3. Confirm EM, PCB, Motor, Stream Function, and IH blocks are masked and use
   the tested application runner.
4. Load and update `matlab/radia_ih_sample.slx`; confirm it has no Python or
   S-function dependency.
5. Run release/deploy gates with `release-qud` when requested.
6. Verify Cubit's toolbar separately.

Do not add PySide6 to normal Radia Python. Do not require MEX unless its
application-specific production gates have passed.
