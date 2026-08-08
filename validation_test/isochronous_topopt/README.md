# Isochronous topology optimization -- validation lane

Executable numerical truth for `radia.isochronous_topopt` (the density-adjoint
HDiv-MMM topology-optimization application; design record:
`docs/hdiv_vim/TOPOLOGY_OPTIMIZATION_ISOCHRONOUS.md`, showcase notebook:
`docs/hdiv_vim/isochronous_topopt.ipynb`).  The fast structural locks live in
`tests/test_isochronous_topopt.py`; this lane re-runs the three stage gates at
the RESEARCH-MESH configurations with golden bands set from the measured
2026-07-28 values.

| Lane | Configuration | Golden band (measured value) |
|---|---|---|
| `test_adjoint_gate_lane.py` | unit ball, maxh 0.35 (~270 tets), log-uniform `s` in [1e-2, 1] | directional adjoint-vs-FD rel < 1e-7 (8.1e-10); dipole-reciprocity vs the independent C++ charge evaluator rel < 1e-7 (1.1e-10) |
| `test_design_loop_lane.py` | sector pole, maxh 0.02 (~194 tets), arc-orbit radial dB_z/dr objective, two mean-B_z arc constraints, 50 % volume budget | J strictly monotone with total gain in [8 %, 30 %] (+16.1 %); violations <= 1.25 band at every accepted iterate (peak 1.06 band); matched-0/1 ersatz band abs < 3 % (+0.69 %); gray-design gap detected (< -50 %, measured -96 %) |
| `test_shape_regen_lane.py` | same sector design -> `radia.topopt_cad` nodal level set -> grid iso STL (watertight, Taubin drift < 2 %) + Cubit `create tri iso` (fixed-caps + free-iso blocks UNION watertight) -> `cubit_stl_to_vol` hex/tet (closure + zero-inverted + boundary-face gates) -> objective re-evaluated on the regenerated Sculpt all-hex body vs the staircase iron-only value | smooth-vs-staircase delta same sign, |delta| < 50 %, hex DemagFactor(z) > 0.2 (measured 0.477 -- locks the 2026-08-08 boundary-face incident where a bare-Sculpt `.vol` with `surfaceelements=0` made the demag solve silently demag-free, 51x off); needs a Cubit license (headless batch). Tet remains an export gate: the measured 16,102-tet solve capped 20,000 CG iterations for mass-Riesz/Jacobi/cluster-tree even at `tol=1e-6`, so its solver scaling is owned by the preconditioner lane. |

Run:

```bash
python -m pytest validation_test/isochronous_topopt -v
```

Wall time ~1 minute on a development host (the solves are warm-started CG on
the build-once charge-Gram H-matrix).  Timings recorded inside the histories
are informal; benchmark-grade study-scale timings run on the quiet compute
hosts per the repository benchmark policy.

Notes fixed by this lane (do not re-walk):

* Constraint LP rows must be normalized to O(1): Tesla-scale rows sit below
  HiGHS's ABSOLUTE feasibility tolerance and read as noise.
* The exact-void verification compares MATCHED 0/1 representations; the
  separate continuous-vs-threshold gap is the gray-design diagnostic (drive
  it down with `penalty`/projection continuation before manufacturing
  conclusions).
* `iron_only_mesh` carries the netgen right-hand-OUTWARD boundary-triangle
  ordering; the opposite handedness flips every surface charge into runaway
  magnetization.
