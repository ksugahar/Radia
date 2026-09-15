# Native evidence refresh, 2026-09-15

**Status: all three reproduced failures resolved; 87/87 MATLAB tests and
46/46 integration/source-freshness checks passed.** The initial investigation
below is historical; the accepted forward run is recorded at the end.

Scope: close the stale motor-angle MATLAB evidence identified during the MCP
package-test boundary review. This is not a full solver release or a live MCP
deployment.

## Reproduction identity

- Native build: clean commit `4cd0a7f4b0a1ca875a3c97a62bf95a8a74656aea`,
  `Build.ps1 -MatlabMexOnly -RequireNativeProvenance` on LAB, exit 0.
- MEX SHA-256: `a9eac0124bc989a77bdfb12a15209bb91c244369bb0f97c2e341f910b303a136`.
- MATLAB test sources: `19f991add` (only per-suite scratch isolation after the
  build commit; numerical source and MATLAB setup/generator unchanged).
- Source archive SHA-256:
  `8d85c13c44fb09159e4f46f46edd0c45a3b8f42134872e46f62e14da90e0e205`.
- Build environment: Python 3.12.10, NGSolve/Netgen 6.2.2606, MKL 2026.1.0,
  pybind11 3.0.2. No editable installation was changed.

## Execution boundary

Dedicated Python MATLAB Engine sessions use `-nodesktop -singleCompThread`.
The foreground supervisor bounds startup at 120 seconds and terminates only
its owned worker tree on timeout. Running desktop sessions are not attached.
Host identity is obtained from the OS; the selected scratch and provenance
environment variables are explicitly passed into the new Engine session.

Both hibino and mdx2 timed out before Engine readiness. Their empty worker logs
were recovered and hash-verified (SHA-256
`e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`).
This is startup failure, not evidence of a numerical failure or a diagnosed
licensing cause. Credentials and runner services were not changed.

LAB Engine startup succeeded. Initial launcher issues (Engine requires
`io.StringIO`, and `COMPUTERNAME` was not inherited) were corrected before the
actual tests. Test results and final disposition are recorded below.

## Result

LAB completed all 87 tests: **84 passed, 3 failed, 0 incomplete**. The generator
correctly returned failure. The unmodified generated result is retained in Git
at commit `e79a6d8e7`, path
`artifacts/native_refresh_20260915/native_motor_angle_family.json`, SHA-256
`576ca559d7c3e768bdfcc43abfcfcda7e046173f23cb84e4d645f3cfffc6c17f`.
It is not promoted into the historical passing artifact directory.

| Failed test | Observation |
| --- | --- |
| `testHACApKHexSelfBlockDirectionalDerivative` | Dilation identity relative discrepancy up to 1.43e-5; central difference relative norm discrepancy 8.61e-6 against 2e-8. |
| `testHACApKNonAffineHexSelfBlockDirectionalDerivative` | Dilation identity relative discrepancy 2.685e-4 against 3e-11. |
| `testApplicationConfigAndFailureArtifacts` | Expected `Coil script` message absent. The source-only staging tree lacks `_radia_pybind`; the Python launch boundary needs a separately validated native Python runtime. This is not proof of a production application failure. |

The first two require investigation of the primal/shape-derivative quadrature
routes in `src/core/rad_hacapk_hdiv.cpp`. Source hashes in the build worktree
and reviewed checkout match. No numerical tolerance or kernel was modified to
turn these failures green. The third is a validation-environment limitation;
do not replace its error assertion with a generic failure check.

At that stage the source-freshness gate remained open; historical native evidence and its
production proof were not restamped. The new proof-generator regression checks
separate version/timestamp fields and the actual linked native JSON hash in a
temporary output, without overwriting that historical proof.

Raw build/Engine logs and transfer inputs are retained on LAB at
`C:/temp/mcp-native-evidence-logs-20260915`, owned by this maintenance task for
the numerical follow-up. Cleanup trigger: accepted replacement evidence and
its committed diagnosis. Dedicated Engine workers have exited; existing
desktop MATLAB sessions were not touched.

Both owned remote staging directories were removed after recovery and evidence
commit; absence was verified on hibino and mdx2. No other job directory was
removed. The clean LAB build worktree
`C:/temp/mcp-native-build-20260915-3c6197863` has the same task owner and cleanup
trigger as the diagnostic logs. A direct isolated `import radia` confirmed the
missing `_radia_pybind` runtime; this must be supplied by a genuine build, not
by weakening the application assertion.

Initial repository integration: 45 passed. Native source-freshness acceptance
was held open until a genuinely passing replacement was generated.

## Accepted forward run

- Native source commit: `1b1687af38ed69d25e12964804a5788e1155f6b6`, clean.
- MEX SHA-256: `bd46e9f152f8801051214616a4bd14819285114ace93074fa48b0727277d8ca6`.
- Python native SHA-256: `7caf525ece2f2b91e1c06cd897b82a800808e54be5f6f930282213c61f517237`.
- Generator commit: `b6406c1cf`; only evidence hashing changed after the native
  build. Kernel and MATLAB tests are unchanged from the native build commit.
- Both `-MatlabMexOnly` and `-RadiaOnly` builds passed. A fresh Python process
  imported Radia 5.0.0 from the selected build tree, not a global editable tree.
- Dedicated LAB R2026a Engine: **87 passed, 0 failed, 0 incomplete**. No tolerance
  was relaxed. Added a non-affine interior-node velocity finite-difference check.
- Python regression: `test_native_production_hex_volume_self_block_derivative_invariants`
  passed, covering translation, scaling and exact symmetry through pybind11.
- BDM1 volume/face self derivatives now differentiate the same pair-domain
  Duffy integral as the primal. Reference Piola charge measures are fixed;
  the derivative kernel is `-dot(X_T-X_S, V_T-V_S)/|X_T-X_S|^3`.
- The passing native JSON is promoted to `annular_motor_dual_lane_v1`, SHA-256
  `88eec885cd8a0f4ddd2bac4f3c79db04025868690b9197eb256bd927b3d5deab`.
  Its production proof is regenerated with separate historical motor/native
  versions. The failed JSON is superseded (Git preserves the investigation).
  Kernel and executed-test hashes now guard evidence freshness, not just the
  MEX gateway source. Integration and source-freshness checks: **46 passed**.
- Current raw logs: `C:/temp/mcp-native-forward-20260915`, owned by this task
  as build/debug evidence until the next native release accepts this fix.
  The current build tree has the same owner and retention condition; it is
  not an editable installation or a rollback distribution.

This closes the three reproduced failures, not a full Radia release, a live
MCP deployment, or the unrelated remote Engine startup investigation.

The obsolete LAB staging directory `C:/temp/mcp-native-evidence-logs-20260915`
could not be removed: the execution policy rejected the recursive cleanup
command before it ran. It remains owned by this task pending permitted cleanup,
not as a rollback target. Current MEX/Python build logs were separately recovered
to `C:/temp/mcp-native-forward-20260915` and their SHA-256 hashes verified.
