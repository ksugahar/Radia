# Multipole-Moment MMM Memory

Date policy: 2026_06_26

This file records short "do not repeat this mistake" notes for the
multipole-moment MMM implementation.  Temporary JSON files are scratch logs;
the durable decisions and key numbers must be written here.

## Engineering Benchmark Range

Decision, 2026_06_26:

- Do not let `mu_r >= 10000` linear high-permeability cases dominate
  engineering conclusions.  They are useful as numerical conditioning stress
  tests, but they are not the main production target for real soft-iron
  devices.
- For performance claims, prioritize nonlinear BH-curve cube/C-yoke runs and
  linear reference cases around ordinary engineering permeability
  (`mu_r ~ 100` to `5000`, with `mu_r ~ 1000` as the common lock/reference
  scale).
- If a preconditioner only helps the `mu_r >= 10000` stress test but slows or
  does not improve the nonlinear/`mu_r ~ 1000` cases, remove it or keep it as
  a separately justified new research branch.  Do not present it as the
  default path.

Two-stage smoke record, 2026_06_26:

| Stage | Case | Method 2 block-Jacobi inner iters | Removed 3-mode coarse correction inner iters | Result |
|---|---:|---:|---:|---|
| A, engineering | linear `mu_r=1000`, 3x3x2 block | 19 | 28 | worse |
| A, engineering | nonlinear BH, 3x4x2 block | 532 | 568 | worse |
| B, stress | linear `mu_r=10000`, 3x3x2 block | 31 | 36 | worse |
| B, stress | linear `mu_r=100000`, 3x3x2 block | 59 | 69 | worse |

Conclusion:

- The 3-mode global dipole coarse correction did not reduce BiCGSTAB
  iterations in either Stage A or Stage B.
- It was deleted from the implementation and public configuration surface.
- The raw scratch JSON was `C:\temp\radia_moment_two_stage_eval_2026_06_26.json`;
  this markdown table is the durable record.

## Method 1 Pure-Hex BiCGSTAB Must Stay Matrix-Free

Decision, 2026_06_26:

- Pure 6-DOF hex multipole-moment method 1 uses `MomentSystemBlock6x6` in split
  form `A(chi)x = Lx + chi*Kx` and stores only the RHS, current iterate, work
  vectors, and 6x6 element-block Jacobi inverses.
- Do not rebuild the dense `SystemMatrix` for pure-hex BiCGSTAB just because
  method 0 still needs dense LU.  That repeats the old comparison-path shortcut
  and hides the real method-1 memory behavior.
- Wedge/pyramid and mixed hex-wedge method 1 also stay matrix-free through
  `MomentSystemBlockAny`; do not bring back a dense variable-DOF method-1 path.
- There is no scalar or identity preconditioner substitute: if an element block
  inverse cannot be built, fail loud and fix the block/preconditioner issue.
- Smoke record: `C:\temp\radia_moment_method1_matrixfree_smoke_2026_06_26.json`
  compares sequential method 0 vs pure-hex method 1 at 24 and 108 DOF;
  external-B relative differences were about `6.3e-13` and `3.4e-12`.
- Follow-up smoke record: `C:\temp\radia_moment_accel_followup_2026_06_26.json`
  compares mixed hex+wedge method 1 against method 0; external-B relative
  difference was about `1.7e-11`.

## Three-Mode Coarse Correction Was Removed

Decision, 2026_06_26:

- The former `SolverConfig(moment_two_level_precond=True)` 3-mode global dipole
  correction was deleted from the implementation and public configuration
  surface.
- The small 108-DOF linear smoke in
  `C:\temp\radia_moment_accel_followup_2026_06_26.json` preserved the solution
  (`relB ~ 1.0e-11`) but did not reduce inner iterations (`14 -> 14`) and
  increased solve time from coarse-space overhead.
- The two-stage smoke summarized above was worse in every small case:
  `19->28`, `532->568`, `31->36`, and `59->69` iterations.
- Passing `moment_two_level_precond` now fails loud.
- Do not reintroduce this 3-mode correction as a casual performance knob.
- A future global/hierarchical preconditioner must be a new design and must
  first pass Stage A engineering cases before any Stage B stress-test gain is
  considered useful.

## Inexact BiCGSTAB Was Removed

Observation, re-confirmed on LAB during the 2026_06_26 acceleration pass:

- Inexact BiCGSTAB can reduce accumulated inner iterations, but the benefit is
  limited compared with nonlinear outer-loop and preconditioner improvements.
- On the 144-DOF nonlinear 3x4x2 lock case, exact method-2 BiCGSTAB used 446
  accumulated inner iterations and 60 outer Picard iterations.
- The tightened inexact schedule reduced this to 344 inner iterations but still
  60 outer iterations, with external-B difference about `4.2e-12`.
- A too-loose early schedule (`1e-4` to `1e-3` class when `bicgstab_tol=1e-9`)
  caused the nonlinear method-2 lock test to miss convergence. Do not repeat
  that schedule.
- Therefore, the implementation was removed from `SolverConfig` and the
  method-2 solve path.  Do not reintroduce it as a default or a casual
  performance knob.

Practical rule:

- Treat inexact BiCGSTAB as a failed optimization branch unless a larger
  benchmark proves that outer convergence and final fields are unchanged.
- Prefer work on nonlinear acceleration, GMRES comparison, and redesigned
  global/hierarchical preconditioning when seeking larger gains.

Relevant measurement:

- `C:\temp\radia_moment_krylov_accel_smoke_2026_06_26.json`
