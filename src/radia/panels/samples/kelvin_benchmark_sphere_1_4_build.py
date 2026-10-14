"""Build Cubit 1/4 sector mesh of magnetic sphere + air + Kelvin and
export Netgen .vol at orders 1, 2, 3.

1/4 reduction: x >= 0, y >= 0 (full z).
  reduction = {"x": "bn=0", "y": "bn=0"}
  offset_dir = "z" (free axis, perpendicular to both reduction planes)

This is the production-validated Kelvin Benchmark sample (verified
2026-04-25): `kelvin_benchmark_sphere_1_4_p2.vol` matches analytical
Hz_origin to +0.71% at p=2 and +0.54% at p=3.  See README at
`validation_test/cubit/kelvin_1_4_p_convergence/`.

Two non-obvious Cubit fixes embedded in `add_kelvin_cubit` and the
shared core (`kelvin_benchmark_sphere_build.py`):

  1. `subtract A from B keep` does NOT carve B in Cubit 2025.3.
     Workaround: subtract WITHOUT keep, then re-create the subtrahend.

  2. Cubit-meshed FaceDescriptor for the `sphere` sideset has
     (domin=air, domout=mag) orientation; the solver therefore uses
     `-specialcf.normal` for the Neumann correction (already the
     OCC convention, so no formulation difference).

Run:
  python kelvin_benchmark_sphere_1_4_build.py [--out-dir DIR]
                                              [--orders 1,2,3]
                                              [--mesh-size 0.025]

Thin wrapper around `kelvin_benchmark_sphere_build.main()` with
`--frac 1_4`.
"""

import sys

from kelvin_benchmark_sphere_build import main


if __name__ == "__main__":
    main(["--frac", "1_4"] + sys.argv[1:])
