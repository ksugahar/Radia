"""Independent high-precision Pade oracle; no native solver or MATLAB needed.

The oracle solves a Taylor-moment Pade system with mpmath, independently of
the production orthogonal projection. Both use the same finite-modal input;
this validates numerical stability, not convergence of the Bessel-mode tail.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import platform
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

import mpmath as mp
import numpy as np
import scipy
from scipy.special import jn_zeros

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out", type=Path, default=Path(__file__).parent / "results/cln_pade_reference.json"
    )
    args = parser.parse_args()
    started = time.perf_counter()
    source = ROOT / "src/radia/maglev/mixed_galerkin/references.py"
    spec = importlib.util.spec_from_file_location("cln_reference", source)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load the analytic reference module: {source}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    rows = []
    # Fix the floating-point modal input, then solve its moments at 100 digits.
    with mp.workdps(100):
        zeros = [mp.mpf(float(x)) for x in jn_zeros(0, 200)]
        moments = [mp.mpf(1)] + [
            (-1) ** k * mp.fsum(4 / x ** (2 * k + 2) for x in zeros) for k in range(1, 21)
        ]
        for order in (1, 4, 10):
            for kind in ("L", "R"):
                p, q = mp.pade(moments, order - (kind == "L"), order)
                for frequency in (1e-3, 1.0, 1e3, 1e6):
                    u = mp.mpc(0, frequency)
                    expected = complex(
                        mp.polyval(list(reversed(p)), u) / mp.polyval(list(reversed(q)), u)
                    )
                    actual = complex(
                        module.Y_cln_pade(
                            complex(u), order, 1.0, 1.0, 1.0, kind=kind, n_modes=200, n_taylor=20
                        )
                        / math.pi
                    )
                    error = abs(actual - expected) / abs(expected)
                    rows.append(
                        {
                            "N": order,
                            "kind": kind,
                            "u_imag": frequency,
                            "reference_ratio": [expected.real, expected.imag],
                            "actual_ratio": [actual.real, actual.imag],
                            "relative_complex_error": error,
                            "passed": error < 1e-8,
                        }
                    )
    git = subprocess.run(
        ["git", "-c", f"safe.directory={ROOT.as_posix()}", "rev-parse", "HEAD"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    import hashlib

    result = {
        "schema": "radia.cln-pade-reference.v1",
        "passed": all(row["passed"] for row in rows),
        "utc": datetime.now(timezone.utc).isoformat(),
        "host": platform.node(),
        "source_base_commit": git.stdout.strip(),
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "elapsed_seconds": time.perf_counter() - started,
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "mpmath": mp.__version__,
        "precision_digits": 100,
        "n_modes": 200,
        "n_taylor": 20,
        "oracle": "100-digit Taylor-moment Pade; shared float64 Bessel zeros only",
        "scope": "Finite-modal Pade evaluation, not infinite-mode or solver validation",
        "relative_tolerance": 1e-8,
        "cases": rows,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(
        f"{sum(row['passed'] for row in rows)}/{len(rows)} passed; max error={max(row['relative_complex_error'] for row in rows):.3g}"
    )
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
