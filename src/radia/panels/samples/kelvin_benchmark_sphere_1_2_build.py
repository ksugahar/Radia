"""Build Cubit 1/2 sector mesh of magnetic sphere + air + Kelvin and
export Netgen .vol at orders 1, 2, 3.

1/2 reduction: y >= 0 (single symmetry plane y=0, BC bn=0 since the
benchmark uniform field is along z and B is parallel to y=0 plane).
Kelvin sphere offset along z (free axis, perpendicular to y=0).

Run:
  python kelvin_benchmark_sphere_1_2_build.py [--out-dir DIR]
                                              [--orders 1,2,3]
                                              [--mesh-size 0.025]

This is a thin wrapper around `kelvin_benchmark_sphere_build.main()`
with `--frac 1_2` -- see that file for the parameterised builder.
"""

import sys

from kelvin_benchmark_sphere_build import main


if __name__ == "__main__":
    main(["--frac", "1_2"] + sys.argv[1:])
