# Induction heating: field, loss and coupling

Radia can resolve conductor skin effects, calculate electromagnetic loss and
couple eddy-current heating to thermal evolution. Look for checked spatial
loss/temperature artifacts, not merely a successful solver exit. The Radia MCP
IH tool family owns the live workflow, supported configurations, dependencies
and validation routes. Simulink's IH blocks are the formal human UI and require
MathWorks' official MATLAB MCP Server.

Start with the result-bearing Python notebooks:

- [ESIM method and benchmark](../ih_esim_benchmark/esim_showcase.ipynb):
  equations, implementation assumptions and saved comparisons.
- [Loop degrees of freedom and cut selection](../cohomology/loop_dof_cut_selection.ipynb):
  why conductor topology matters to current representation.
- [Dowell and surface-impedance comparison](../peec/dowell_surface_impedance_demo.ipynb):
  analytical loss models and their frequency/geometry limits.

These are discovery and reproducibility pages, not notebook workbenches.
Python implementations are LLM-driven through MCP. Electromagnetic examples
alone do not certify a coupled transient thermal application. Numerical
acceptance remains in [IH validation](../../validation_test/induction_heating/).

The former source/hash catalogs are preserved in
[documentation maintenance](../../validation_test/documentation_maintenance/README.md);
they are not public numerical demonstrations.
