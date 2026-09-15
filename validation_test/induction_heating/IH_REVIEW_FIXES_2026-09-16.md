# IH review corrections (2026-09-16)

Scope: review of main 6793633, implemented on 7c931c510 plus this change.
This is partial acceptance, not a release certificate.

## Findings and disposition

- A/B: confirmed. Thermal incorrectly rotated a material-frame state after
  Eddy's source/material mapping. Removed that extra rotation. A dependency-free
  native regression failed before the fix and passes after it. Integer and
  half-cell full turns heat all eight cells from 300 to 301 K; passive material
  temperatures remain fixed. The same full-turn oracle passes through MEX.
- C: linear interpolation remains explicitly diffusive; conservation is only
  of a weighted integral. No claim of peak preservation or high-order rotation.
- D: reject ill-conditioned weighted-integral corrections using a relative
  cancellation threshold. The signed-field reproducer now throws; nonfinite
  inputs/integral overflow also fail.
- E/F/G: document the diagonal loss representation, frozen unit-current loss
  generator, unavailable production nonlinear BH generator, and dense slope
  path cost. No unvalidated tolerance-based freezing of EM operators was added.
- H: monitor caches its validated config at Start rather than validating all
  CSR operators each Outputs call. FE temperature-statistics evaluation still
  validates its own evaluation operator; further optimization is not claimed.
- I: physical nonfinite signals fail. Review's interpretation of output 4 was
  incorrect: addIHMonitorBus wires it to header.error_code, not a sample count.
  Preserve zero on success and correct the comment. Output 5 already counts
  samples. NaN geometry revision is an intentional cache-invalid marker, not
  a physical signal failure.
- J: Eddy reset added at native, MEX and S-Function boundaries.
- K/L: document weights as representation validation (M already owns weights);
  compare full mass/stiffness row pointers and reject nonzero initial CSR offset.
- M: warm start CG from accepted temperature, with RHS-relative stopping bound.
  No preconditioner/performance claim added.
- N: current flattened ABI remains explicitly row-major. Document A(:) hazard;
  mandatory layout-tag migration is NOT implemented in this change.
- O: cleanup exceptions emit warnings, and mexLock is outside registry mutex.
- P: retain C:/temp per Windows machine policy; portable test-root migration
  is not included.

## Checks

- MSVC dependency-free native regression: pass; CTest target radia.ih_rotation.cpp.
- MATLAB R2026a Update 3, dedicated Python Engine, rebuilt MEX: 14 passed,
  0 failed, 0 incomplete across native S-Function integration and model tests.
  The initial tracked-model assertion failed while querying ObjectParameters.
  Independent official Toolkit inspection showed all three controls were present
  in the saved model. Test the mask contract directly with MaskNames and assert
  the exact resolved FileName. No missing SLX control was established; the earlier
  missing-control diagnosis was incorrect. Generated and tracked model tests pass.
- Application interface manifest: 5 passed. Application block audit: CLEAN.
- test_simulink_application collection blocked: isolated checkout lacks the
  Python _radia_pybind binary. No other checkout's binary was substituted.
- Standard MEX build compiled/linked, but post-build copy failed because the
  existing destination was loaded. Tests used the new binary in an isolated
  directory; no borrowed MATLAB session was stopped or unloaded.
- Official MATLAB MCP failed to attach, but the installed official Toolkit
  executed through the dedicated Python Engine: library prerequisite passed,
  model_read succeeded and model_check reported healthy. No SLX edit was needed;
  no ZIP/XML repairs were performed. Exact-package cross-host gates, full visual
  release acceptance and publication remain open.

See ih_review_fixes_20260916.json for binary identity and acceptance totals.
