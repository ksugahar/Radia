"""Build Cubit 1/8 sector mesh of magnetic sphere + air + Kelvin and
export Netgen .vol at orders 1, 2, 3.

1/8 reduction: x >= 0, y >= 0, z >= 0.
  reduction = {"x": "bn=0", "y": "bn=0", "z": "ht=0"}
  offset_dir = "x" (no free axis exists; the kelvin_far cut becomes
  the Dirichlet "infinity" boundary)

ACADEMIC-ONLY -- DO NOT SHIP THE .vol IN panels/samples/.

The 1/8 case is a sphere-symmetric academic benchmark.  Realistic
electromagnet geometries (C-magnets, dipoles, quadrupoles) at most
have 1/2 or 1/4 symmetry; 1/8 only applies to spheres / cubes /
fully-symmetric primitives, which are better served by analytical
checks or fixed-source Radia field evaluations than by FEM-Kelvin.

The build script is kept here so future research can pick it up,
but the EM panel itself is fully covered by the 1/2 and 1/4 samples
(`kelvin_benchmark_sphere_1_2.vol`, `kelvin_benchmark_sphere_1_4.vol`).

Solver status (2026-04-26): Hz_origin -> 0 in the continuous limit
(verified by h-refinement); the failure is formulation-level, not
discretization.  rho_min cap and Dirichlet variants do NOT close
the gap.  See `memory/feedback_kelvin_1_8_not_a_real_use_case.md`
and `memory/feedback_kelvin_1_8_blocker.md` for the investigation.

Run:
  python kelvin_benchmark_sphere_1_8_build.py [--out-dir DIR]
                                              [--orders 1,2,3]
                                              [--mesh-size 0.025]

Thin wrapper around `kelvin_benchmark_sphere_build.main()` with
`--frac 1_8`.
"""

import sys

from kelvin_benchmark_sphere_build import main


if __name__ == "__main__":
    main(["--frac", "1_8"] + sys.argv[1:])
