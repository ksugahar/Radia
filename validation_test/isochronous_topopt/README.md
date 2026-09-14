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
| `test_shape_regen_lane.py` | same sector design -> `radia.topopt_cad` nodal level set -> grid iso STL + Cubit-native iso surface -> gated hex/tet `.vol` -> objective re-evaluated on the regenerated all-hex body | watertightness and volume-drift bands; closure, inversion, and complete-boundary gates; smooth-vs-staircase sign and magnitude band; nonzero demagnetizing response. Needs a Cubit license; tet solver scaling remains owned by the preconditioner lane. |
| `test_hex_native_design_lane.py` | the complementary Sculpt + HDiv-MMM half: sector domain STL -> `cubit_stl_to_vol` hex (size 0.012, 360 hexes, min quality 0.540) -> the SLP design loop and the hex path of `iron_only_mesh` / `verify_design_iron_only` | all mesh gates incl. `boundary_faces_ok`; accepted iterates monotone with gain in [0.5 %, 15 %] (+8.2 % in 30 iterates measured 2026-08-10, after the band-edge acceptance-zone fix un-jammed the SLP); iron extraction stays hex; matched-0/1 ersatz bands abs < 20 %. Requires a Cubit license and runs headless (~3-4 min). |
| `test_vfrac_route_lane.py` | the Sculpt-native volume-fraction regeneration + solver chain: sphere body -> `write_vfrac_exodus` -> `cubit_vfrac_to_vol` (standalone `sculpt.exe --input_vfrac`) -> `hdiv_vim_demag_eval` | bridge gates all green with closure < 1e-2 (measured 2e-3; the STL route's gate is 3e-2); air-mesh-free demag factor inside [0.30, 0.37] against the sphere's closed-form 1/3. Requires a Cubit license (~2 min). |
| `test_vfrac_periodic_lane.py` | Sculpt `--periodic` on the volume-fraction route: a simple-cubic array of overlapping spheres (R = 0.55 L, minimum-image sampling) meshed as a periodic RVE, with the same input meshed plainly as the negative control | the writer's `bounds=` spans exactly one period with no guard ring; the periodic cell's node set is invariant under the period (matched +/- counts equal, seam fill rate >= 0.98 -- measured 1.000 on all three axes) while the plain mesh fills only 0.061 of its seam; solver gates (closure, inversion, outer skin) unaffected. Requires a Cubit license (~1 min). |
| `test_vfrac_multimaterial_lane.py` | the same chain with a PARTITIONED two-material body (core r<=0.3 inside shell r<=0.6, union = the single-material sphere): conformal Sculpt interface, named `.vol` regions, per-region `mu_r` | gates green with the OUTER faces matching the topological skin and interface faces > 0 (measured 894 + 190); core/shell share interface nodes; uniform `mu_r` recovers the sphere's 1/3 (measured 0.327) with per-region mean-magnetization ratio ~1 (1.001), while a 200x contrast drives that ratio below 0.5 (0.006) -- the demag factor is geometric and stays 0.32734 in both, so it cannot be the discriminator. Requires a Cubit license (~5 min). |

Run solver-heavy lanes only on an idle compute host under the repository host
policy. Do not run this whole directory on LAB/100 merely because Cubit is
installed there. The former one-minute development-host estimate does not apply
to shape reacceptance with resolved facet boundaries.

For the post-Taubin shape lane, use the staged driver:

```bash
python validation_test/isochronous_topopt/run_shape_reacceptance.py prepare --directory C:/temp/shape-run
# Transfer prepare.json, prepare.state.json, both STLs and design_lsd.exo to LAB.
python validation_test/isochronous_topopt/run_shape_reacceptance.py mesh --directory C:/temp/shape-run --command-plugin-directory "C:/Program Files/Coreform Cubit 2025.12/bin/plugins" --tet-curve-interval 16
# Transfer the resulting directory back to the same isolated compute environment.
python validation_test/isochronous_topopt/run_shape_reacceptance.py evaluate --directory C:/temp/shape-run
# After recovering the completed directory, independently audit receipt consistency.
python validation_test/isochronous_topopt/run_shape_reacceptance.py audit --directory C:/temp/shape-run
```

Use a fresh owned directory; a repeated or failed phase is deliberately refused.
For Python-only CAD experiments, all compute phases use the same explicit
`--topopt-cad` source, and the record identifies that overlay separately from
the installed numerical wheel. It is not a release acceptance of a rebuilt wheel.
The prepare phase records the staircase scalar from the same optimization that
produced the STL; re-running optimization is not a valid way to reconstruct its
identity because small floating-point differences can change decimation.

The interval override changes the original mesh-generation conditions: it resolves
all imported facet boundary curves, including tiny sliver surfaces that the old
global size missed. No CAD face is deleted, relabelled away, or exempted from
`check-vol`. All commands and clean process exits are recorded. The original
field tolerances remain unchanged. A successful mesh receipt alone is not field
acceptance: require `field.json` with `completed=true`, its matching completed
evaluation state, and the verified prepare/mesh parent receipt chain.

### Post-Taubin acceptance, 2026-09-14

`shape_reacceptance_20260914.json` records a completed mdx2 run against candidate
wheel source `59b094d8ed2c19e35967fd7631f3a1473eb3d200`, with the explicit
`topopt_cad.py` overlay identified inside the result. This is acceptance of this
shape lane, not acceptance of every solver or of a new binary release.
The native/NGSolve objective discrepancy was at most `6.272e-10`; independent
direct-field reciprocity discrepancy was at most `4.038e-8`. The refined hex
objective differed from the tetrahedral reference by 2.2805%, improving on the
coarse hex's 16.3797%; it differed from the original staircase by -2.9249%.
All pre-existing field tolerances passed without relaxation.

The first-party fixture `fixtures/shape_reacceptance_20260914.zip` contains the
exact generated geometry, solver meshes, mesh reports, original field result, parent receipts and
original source variants. No externally downloaded dataset is included. To
audit it, extract into a fresh directory and run the `audit` phase. Use the
archived `field.json`: Git may normalize line endings in the readable tracked
JSON copy, whereas the receipt digest identifies original bytes. The separate
`shape_reacceptance_20260914_audit.json` records post-execution validation;
`numerical_recomputed=false` means it does not impersonate a new solver run.
Original numerical driver SHA is preserved rather than replaced by the later
validator's SHA. Full console logs and failed mesh trials are retained on LAB,
not published as repository log files.

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
