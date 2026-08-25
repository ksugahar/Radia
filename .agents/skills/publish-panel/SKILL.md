---
name: publish-panel
description: Publish Radia's Simulink application library after verifying its five application blocks and tracked MATLAB-only samples. Use simulink-app-health before deploy/release-quad; Cubit's toolbar remains a separate deployment surface.
---

# Publish Application Blocks

1. Run `simulink-app-health`.
2. Regenerate and inspect `matlab/radia_simulink_library.slx`.
3. Confirm EM, PCB, Motor, and Stream Function use the tested application
   runner. Confirm IH contains the native Eddy/Thermal MEX S-Functions.
4. Load, update, and simulate `matlab/radia_ih.slx`; confirm its direct
   Thermal-to-Eddy feedback and native-only backend. Treat the first package as
   a preassembled-operator preview until native `.vol` assembly lands.
5. Build the versioned archive and run `release_quad simulink-candidate` for all
   four machines. `release_quad done --simulink-package <zip>` must pass before
   GitHub publication.
6. Run the remaining release/deploy gates with `release-quad` when requested.
7. Verify Cubit's toolbar separately.

Do not add PySide6 to normal Radia Python. IH is the explicit native-MEX
exception and must pass its application-specific production gates.
