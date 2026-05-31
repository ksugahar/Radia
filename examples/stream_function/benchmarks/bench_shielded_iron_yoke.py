"""STUB: Shielded coil with iron back plate via Radia MMM kernel (TODO).

This stub is the material-kernel extension demo.  The (A) callback
contract lets us swap the free-space Biot-Savart kernel for Radia's
MMM/MSC kernel that handles magnetic materials (iron yoke, permanent
magnet shield, conducting shield via SIBC, ...).  Implementing this
exercises the kernel-agnostic interface end-to-end and produces
numbers that can be reported alongside the free-space benchmarks.

Next steps:
  1. Define source plane + iron back plate geometry (e.g., a 5 mm steel
     plate behind the source plane).
  2. Build Radia container with the iron yoke as ObjRecMag with
     MatSatIsoTab BH curve.
  3. Replace bench's entry(i, j) with one that calls rad.Solve() then
     rad.Fld() for each (target, basis) pair.  Expensive per entry --
     ACA+ becomes useful at this scale.
  4. Run our pipeline, log JSON.

Expected behaviour: per-entry ~100x slower than free-space, ACA+ saves
roughly the M / k_aca ratio (~25x at our M, growing with M).  The iron
back plate should mirror the source's field toward the target (like a
magnetic ground plane), improving uniformity at the same coil power.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from benchmark_framework import Benchmark, BenchmarkSpec


class ShieldedIronYokeBenchmark(Benchmark):
    @property
    def spec(self):
        return BenchmarkSpec(
            name="shielded_iron_yoke",
            reference="Material-kernel demo: Radia MMM via the (A) callback",
            target_kind="uniform_Bz_with_iron",
        )

    @property
    def baseline(self):
        # Compare to free-space (no iron) baseline running through the
        # same pipeline -- shows the iron back plate improvement.
        return None

    def solve(self):
        raise NotImplementedError(
            "Shielded iron yoke benchmark not implemented; "
            "see ../README.md TODO.")


if __name__ == "__main__":
    raise SystemExit("Shielded coil stub; see ../README.md TODO.")
