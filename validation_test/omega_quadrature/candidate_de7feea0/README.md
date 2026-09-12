# ESRF 6 candidate-wheel diagnostics: HOLD

These are preserved research observations, not three-engine field acceptance,
not a published-wheel certificate, and not evidence that all HDiv-MMM cases are
complete. This is linear iron with mu_r=1000 on the public ESRF example-derived
CoilBuilder model, not the original example's nonlinear material law. No CST,
TOSCA, or JMAG result, commercial project file, license key, or solver binary is
included. The JSON contains laboratory paths and host/runtime provenance.

## Reproduction Identity

- Candidate source: `de7feea0dea67f0a15209c1b64a39f53fec47adb`.
- Unpublished wheel: `radia-4.95.91-cp312-cp312-win_amd64.whl`.
- Wheel SHA-256: `4952ba92cc0502e81949c35745108f6da499d6c6eb5bf1055a46de1be3a9ad4a`.
- Native SHA-256: `5d8fe9320b6776ed5712ca7158b9878994a827d99d5725278263734dbb8a5ac8`.
- Original diagnostic lane: `6d3af08f3`; original file hashes are in preflight/result JSON.
- Mesh SHA-256: `dfc12b84f80db1fa17fb5012b6f87c072e8aa1fb3c085b61269f047a879a54ac`.
- Mesh: 134686 curved order-2 tetrahedra. `check_vol.json` passed the strict
  sibling `../esrf6_labels.json` contract. The mesh itself is not included.

`threads8/telemetry.json` records the exact command. It uses the isolated
hibino candidate venv, wheel mode, p=1, assembly bonuses 4/8, evaluation order
16, fixed source order 3 and eight threads. Use the pinned lane to reproduce
the historical result; running the current gate is a separate reanalysis.
No new computation was launched while preserving these files.

## Completed Eight-Thread Run

The process finished with exit 2 after 20115.007 s. Both solves completed;
bonus 4 failed the relative phi_total residual gate. Bonus 8 passed the
recorded diagnostic gates. This does not establish assembly/evaluation
quadrature convergence, p convergence, or three-formulation field agreement.

| Assembly bonus | Solve time (s) | Solve plus audit (s) | Diagnostic gates |
| --- | ---: | ---: | --- |
| 4 | 211.534 | 9523.312 | FAIL: relative phi_total residual |
| 8 | 706.235 | 10513.505 | PASS |

At bonus 4, phi_total has residual norm about 7.65e-17 and RHS norm about
2.33e-15, giving relative residual 0.0328. This is retained as a failed
normalization gate, not silently converted to a pass. A justified scaled
residual criterion remains separate work. The energy-identity gates passed.

The OS peak working set and periodically sampled process-tree RSS in telemetry
are different measurements. They must not be substituted for one another.
Most elapsed time here is auditing, not the linear solve.

## Initial Failure and Heap Probes

The 60-thread attempt failed after 287.432 s during independent quadrature
with `Local Heap overflow` (`integrate-lh`, 166666 bytes). No result JSON was
produced. Its telemetry is preserved; the raw temporary log is not published.
`heap_probe.json` and `heap_probe_threads.json` document the small reproduction:
changing SetHeapSize did not change that overload's local heap, whereas one
and eight threads completed. This observation is specific to NGSolve 6.2.2606
and the tested integration overload, not a general thread-safety conclusion.

## Preservation Boundary

Eight JSON files, including the label contract, were copied unchanged from the
local diagnostic staging tree and SHA-256 compared before Git normalization.
The source result JSON SHA-256 was
`db05a5c6fd9a43babe03768c7e553fac13d630bae67edbcbeac7cf4a9ec2b5fd`.
Existing Python/native hashes and HOLD labels were
not rewritten to describe the current checkout. The current tests recheck
identity consistency and the historical gate decisions without solving.
