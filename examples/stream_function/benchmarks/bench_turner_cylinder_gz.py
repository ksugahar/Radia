"""STUB: Benchmark vs Turner cylindrical Gz gradient coil (TODO).

Reference: Turner, R., "A target field approach to optimal coil design",
J. Phys. D: Appl. Phys. 19, L147 (1986).  IEEE TMI 5 (1986) follow-up.

Turner gives ANALYTICAL SFD expressions for a cylindrical source surface
with a Gz gradient target inside the cylinder.  Our existing
demo_sf_to_peec_gz.py already does cylindrical Gz; this benchmark just
loads Turner's published target spec and reproduces.

Next steps:
  1. Extract Turner's analytical SFD (Bessel function expansion).
  2. Compare our discretised SFD to the analytical reference.
  3. Compare wire pattern + DSV uniformity vs Turner's published numbers.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from benchmark_framework import Benchmark, BenchmarkSpec


class TurnerCylinderGzBenchmark(Benchmark):
    @property
    def spec(self):
        return BenchmarkSpec(
            name="turner_cylinder_gz",
            reference="Turner, J. Phys. D 19, L147 (1986); IEEE TMI 5",
            target_kind="gradient_Gz",
        )

    @property
    def baseline(self):
        return None

    def solve(self):
        raise NotImplementedError(
            "Turner cylinder Gz benchmark not implemented; "
            "wraps existing demo_sf_to_peec_gz.py with Turner's "
            "analytical SFD as the baseline.")


if __name__ == "__main__":
    raise SystemExit("Turner benchmark is a stub; see ../README.md TODO.")
