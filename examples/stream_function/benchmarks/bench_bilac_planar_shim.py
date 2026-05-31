"""STUB: Benchmark vs Bilac et al. planar shim coil design (TODO).

Target spec to be extracted from:
    Bilac, V., et al., "Optimal shim coil set for ..." (citation TBD).
    Planar shim coil, B0 + Z gradient + ZZ second-order shim.

This stub fails explicitly with NotImplementedError so the benchmark
suite's CI doesn't report a silent pass.  See ../README.md TODO row.

Next steps:
  1. Locate the Bilac paper PDF (W:\01_paper\?).
  2. Extract target field formula, source plane dimensions, target DSV.
  3. Run our SF design + Path-A + deformation, log JSON result.
  4. Compare published baseline numbers (RMS, wire length, inductance).
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from benchmark_framework import Benchmark, BenchmarkSpec


class BilacPlanarShimBenchmark(Benchmark):
    @property
    def spec(self):
        return BenchmarkSpec(
            name="bilac_planar_shim",
            reference="Bilac et al. (TBD); planar shim coil design",
            target_kind="uniform_Bz",
        )

    @property
    def baseline(self):
        return None      # awaiting literature extraction

    def solve(self):
        raise NotImplementedError(
            "Bilac benchmark not implemented; see ../README.md TODO.")


if __name__ == "__main__":
    raise SystemExit("Bilac benchmark is a stub; see ../README.md TODO.")
