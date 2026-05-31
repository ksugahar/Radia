"""STUB: Benchmark vs Lemdiasov & Ludwig 2005 target field MRI coils.

Reference: Lemdiasov, R.A., and Ludwig, R., "A stream function method
for gradient coil design", Concepts Magn Reson Part B 26B(1), 67-80 (2005).

Detailed target field formulation with explicit published numerics for
multiple coil designs.  Good "next-after-Turner" benchmark for our
machinery.

Next steps:
  1. Download the paper (W:\01_paper\).
  2. Extract published target spec + their reported metrics.
  3. Wrap our pipeline + JSON-out for direct comparison.
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from benchmark_framework import Benchmark, BenchmarkSpec


class LemdiasovLudwig2005Benchmark(Benchmark):
    @property
    def spec(self):
        return BenchmarkSpec(
            name="lemdiasov_ludwig_2005",
            reference="Lemdiasov & Ludwig, Concepts Magn Reson Part B 26B (2005) 67",
            target_kind="gradient_Gz",
        )

    @property
    def baseline(self):
        return None

    def solve(self):
        raise NotImplementedError(
            "Lemdiasov-Ludwig 2005 benchmark not implemented; "
            "see ../README.md TODO.")


if __name__ == "__main__":
    raise SystemExit("Lemdiasov-Ludwig benchmark stub; see ../README.md TODO.")
