"""STUB: Head-to-head comparison vs CoilGen open-source OSS (TODO).

Reference: Schwartz, A., et al., "CoilGen: an open source MATLAB-based
toolbox for the design of coils for magnetic resonance" (GitHub:
github.com/Philipp-MR/CoilGen).

CoilGen is the only other open-source SF-based coil designer.  This
benchmark runs both ours and CoilGen on the same target spec and
compares (RMS, wire length, inductance, compute time).

Next steps:
  1. Install CoilGen + MATLAB / Octave.
  2. Pick an example case CoilGen ships with (e.g., a Gx gradient or
     uniform-Bz shim).
  3. Reproduce the same target in our framework.
  4. Side-by-side metrics table.

Expected outcome: our framework wins on flexibility (regularisation
choices, deformation outer loop, kernel-agnostic for materials) but
CoilGen may be faster / more polished on the standard MRI gradient
case where they have years of tuning.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from benchmark_framework import Benchmark, BenchmarkSpec


class CoilGenHeadToHeadBenchmark(Benchmark):
    @property
    def spec(self):
        return BenchmarkSpec(
            name="coilgen_headtohead",
            reference="Schwartz et al., CoilGen "
                       "(github.com/Philipp-MR/CoilGen)",
            target_kind="gradient_Gz",
        )

    @property
    def baseline(self):
        return None

    def solve(self):
        raise NotImplementedError(
            "CoilGen head-to-head benchmark not implemented; "
            "needs CoilGen+MATLAB install.  See ../README.md TODO.")


if __name__ == "__main__":
    raise SystemExit("CoilGen benchmark stub; see ../README.md TODO.")
