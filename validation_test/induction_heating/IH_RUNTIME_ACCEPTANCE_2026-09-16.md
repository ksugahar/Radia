# IH runtime acceptance follow-up (LAB, 2026-09-16)

This is a scoped acceptance record, not a complete IH release certificate.
The numerical implementation and executed public thermal notebook were merged
through PRs #280 and #281 respectively.

## Cubit export diagnosis: execution route, not an IH solver defect

The earlier persistent `cubit_exec` attempt used Cubit's bundled-Python daemon.
That route accepts native meshing commands but does not load C++ command plugins.
Consequently `export netgen` failed at parsing, before the exporter or solver ran.
The existing `run_headless_journal` implementation explicitly documents this
distinction. The previous attempt selected the wrong execution route; it did
not prove that the deployed exporter was broken.

The same saved cylinder checkpoint was then opened through the existing checked
headless journal runner, using `coreform_cubit.com`. With the deployed plugin
directory explicitly selected and `-noinitfile`, the run exited **0**, loaded
`cubit_mesh_export.ccm`, and exported a curved order-2 mesh. No installation,
user startup file, or solver source was changed. An earlier run that included
user initialization produced the file but exited 2; it is not the accepted run.

| Check | Observed result |
|---|---:|
| Cylinder radius / height | 0.025 m / 0.025 m |
| Volume elements | 181 tetrahedra |
| Geometry order at export | 2 |
| CAD volume relative error | -0.0162604% |
| CAD surface area relative error | -0.0173171% |
| Minimum scaled Jacobian | 0.313920 (approximately) |
| Production `ih_workpiece_v1.json` strict label contract | passed |
| Canonical `check-vol`, including curved mapping quality | passed |

The same coarse mesh exported at order 1 failed the 1% CAD geometry gate
(volume -2.56%, area -2.02%). That negative result was preserved; the threshold
was not relaxed. Curvature was applied by the exporter while CAD was live,
not by calling `Curve()` after loading a `.vol`.

The legacy MCP `cubit_check_vol` call timed out after 300 seconds. The canonical
`python -m cubit_mesh_export.check` CLI completed and produced the accepted
report with the production label contract. This is CLI checker acceptance,
not a claim that the legacy MCP transport was repaired.

The recovered files live under
`Radia/release-quad/ih-p2-candidate-7aec31d72/runtime-20260916` on the lab NAS;
`SHA256SUMS.json` binds the checkpoint, VOL, CAD sidecar, export result and
checker reports. These are operational artifacts, not committed log files.

## MATLAB Engine: environment gate remains open

Both mdx runner services were observed running as `.\\Administrator`, and both
Python 3.12 interpreters imported MATLAB Engine. No shared Engine or running
MATLAB/ServiceHost process was observed during the read-only inventory.

The mdx2 ServiceHost log from the preceding failed run records online licensing
error 5002, `Authentication canceled`, at 2026-09-15 23:21:04. This establishes an
unfinished authentication flow, not who canceled it or why. Authentication must
be completed by the user before another owned startup/calculation/shutdown
probe and the exact-package gate. No credential store or runner was modified.
The mdx1 startup failure is not assigned the same cause without matching evidence.

## Simulink visual observations

Python Engine attached by the explicitly discovered shared name to the existing
LAB MATLAB R2026a Update 3 (PID 26692). It confirmed that `radia_ih` resolved to
`C:/temp/radia-ih-p2-integration/matlab/radia_ih.slx`, with `Dirty=off` and a stopped
simulation. The borrowed MATLAB was not restarted or closed.

The supported Windows Computer Use route then displayed the actual Simulink
window, including title bar, left/right dock strips and status bar, and the
IH Parameters mask. The inspected Radia-authored labels were English and no
replacement glyphs or question-mark corruption were observed. These screenshots
are present in the task record. This observation is separate from numerical
acceptance and does not constitute a fresh installed-package simulation.

## Remaining acceptance

- Complete mdx authentication/startup and repeat the four-host exact-package gate.
- Complete clean installed-package open/run and the remaining UI interactions.
- Complete the unmocked CAD/VOL-to-EM-to-axisymmetric-heat chain. The checked
  cylinder here is the **3D EM workpiece input**, not a 2D thermal mesh.
  A supported 2D `(r,z)` export/fixture route and coil input still need to be
  established. Do not reinterpret this 3D file as a 2D mesh or substitute the
  earlier in-memory thermal test for real-file acceptance.

No new release was published by this follow-up.
