---
name: pyside6-health
description: Guard the boundary between retired Radia desktop panels and Coreform Cubit's private PySide6 runtime. Current Radia application interfaces are Simulink blocks; use simulink-app-health. Never uninstall Cubit's bundled PySide6.
---

# PySide6 Boundary Guard

Radia's production application interface is Simulink. Normal Radia Python on
LAB, 100, mdx, and hibino must not gain a PySide6 dependency.

Coreform Cubit's embedded Python owns its bundled PySide6 and the in-Cubit
export toolbar may use it. Never remove or alter that private runtime as part
of Radia cleanup.

Use `simulink-app-health` for the five application blocks,
`ipynb-gui-health` only for IH's temporary notebook comparison, and
`cubit-plugin-install --verify-only` plus `cubit-smoke-test` for Cubit.
