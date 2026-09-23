# Standalone vacuum source reproduction

This is NOT a three-engine acceptance run or a matched-error performance
comparison. The HDiv material solver excludes relative permeability one;
only the existing mixed Omega adapter is exercised here, without calling
the three-engine driver's main acceptance route.

On hibino, retained Radia 5.0.0 wheel, NGSolve 6.2.2606, eight threads,
checked coarse Kelvin mesh, source and response order two, bonus quadrature
four, relative permeability one:

- Raw vector relative RMS against the prescribed coil field: all points
  0.0120473258%; core abs(x) <= 10 mm: 0.0125539701%.
- Main DOF 41305, adapter runtime including source preparation 676.5436092 s.
- Source preparation 369.7287686 s; 4277500 cache hits, zero misses.
- Sampled process RSS maximum 1018044416 bytes at 50 ms intervals, not an
  isolated allocation peak; sampling can miss peaks.
- Kelvin tangential trace residual 3.9678932%, below the existing 5% gate.
  This gate is not a physical-field error tolerance.

No threshold was selected after observing this result. Small vacuum error
does not certify magnetic-material accuracy, identify the cause of the
remaining magnetic HDiv/Omega difference, or prove split invariance for all
meshes/orders. Do not subtract this error from the magnetic comparison:
the material, mesh level and symmetry processing differ.

The runner exited zero. Raw logs, wheel, checked inputs, driver and cache
wrapper are retained in LAB archive:
S:/Radia/validation_artifacts/hdiv_vacuum_20260924/hdiv-vacuum-20260924-recovery.zip
SHA256 8870d0675c34b8c907d4a832992914a80ddcb2b1410947785842de2a366bd238
(remote and recovered archive matched).

Wheel SHA256 41b5f5ec02c5c26ad47147da7e4496525e3e2066216d85d787487ecc253023e6.
Driver/cache wrappers are unchanged from the full-finer evidence d3e74c962;
their hashes and loaded runtime paths are in environment.json.
Foreground command: ssh hibino 'pwsh -NoProfile -File C:/temp/hdiv-vacuum-20260924/run_vacuum.ps1'.

Next: retain the existing manuscript/slide speed-claim limits. Investigate
compatible source/response refinement with magnetic material or controlled
gap discretization before selecting equal-error performance points.
