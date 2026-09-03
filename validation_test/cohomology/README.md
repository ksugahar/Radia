# Cohomology validation

This directory owns durable numerical evidence for the cohomology algorithms.
Public notebooks under `docs/cohomology/` embed their demonstration results and
WebGUI scenes, but do not write adjacent JSON sidecars.

`loop_dof_cut_selection_data.json` records the accepted genus-one torus cut
selection. Its static test derives the topology and winding checks again from
the stored record so edited evidence cannot silently weaken the gate.
