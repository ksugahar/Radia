---
name: panel-preview
description: Compatibility note for previewing Radia interfaces after the PySide era. Inspect the generated Simulink library and masked application blocks; use ipynb-gui-health only for the temporary IH notebook.
---

# Interface Preview

Build `matlab/radia_simulink_library.slx` in `C:\temp`, open the library, and
inspect the five masked application blocks. The production artifact is the
Simulink library, not a Qt screenshot or docs notebook.
