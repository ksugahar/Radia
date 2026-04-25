"""Build Cubit 1/8 sector mesh of magnetic sphere + air + Kelvin and
export Netgen .vol at orders 1, 2, 3.

1/8 reduction: x >= 0, y >= 0, z >= 0.
  reduction = {"x": "bn=0", "y": "bn=0", "z": "ht=0"}
  offset_dir = "x" (no free axis exists; the kelvin_far cut becomes
  the Dirichlet "infinity" boundary)

WARNING: 1/8 mesh build succeeds at machine precision, but the
Omega-Reduced + Kelvin solve is currently NOT validated -- the
kelvin_far Dirichlet plane passes through r' = 0 (Kelvin centre)
where the reluctivity Mu = mu0*(R/r')^2 is singular, and the
ill-conditioned linear system gives Hz ~ 0 instead of the analytical
3/(mu_r+2) * H0.  Track this in `memory/feedback_kelvin_1_8_blocker.md`.

Until the singular-reluctivity formulation is fixed, prefer 1/4
(`kelvin_benchmark_sphere_1_4_build.py`) which has a free offset axis
and matches analytical to <1% at p=2.

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
