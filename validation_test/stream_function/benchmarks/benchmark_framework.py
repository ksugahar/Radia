"""Common framework for the SF coil design benchmark suite.

Each benchmark subclass defines:
  - ``target_spec``: dict with reference + target field + geometry
  - ``solve()``:    runs our pipeline, returns the designed coil
  - ``metrics()``:  computes RMS, p2p, length, etc.
  - ``baseline``:   dict with the published reference numbers

The framework handles JSON serialisation, comparison, and regression
locking so benchmarks reproduce identically across sessions.
"""
from __future__ import annotations

import abc
import argparse
import json
import sys
import time
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Any

import numpy as np


@dataclass
class BenchmarkSpec:
    """Frozen description of a benchmark's target / reference."""
    name: str
    reference: str         # citation
    target_kind: str       # "uniform_Bz" / "gradient_Gz" / "gradient_Gx"
    target_params: dict[str, Any] = field(default_factory=dict)
    geometry: dict[str, Any] = field(default_factory=dict)


@dataclass
class BenchmarkResult:
    """One run's measured metrics + baseline + ratios."""
    name: str
    reference: str
    target_spec: dict
    our_result: dict
    published_baseline: dict | None
    ratio: dict | None
    elapsed_s: float


class Benchmark(abc.ABC):
    """Abstract base class for a single benchmark case."""

    @property
    @abc.abstractmethod
    def spec(self) -> BenchmarkSpec:
        ...

    @property
    @abc.abstractmethod
    def baseline(self) -> dict | None:
        """Published numbers from the reference, or None if not extracted."""
        ...

    @abc.abstractmethod
    def solve(self) -> dict:
        """Run our pipeline; return result dict with at least 'rms',
        'p2p_mean', 'wire_length_m', 'n_contours' keys."""
        ...

    def run(self, json_out: Path | None = None) -> BenchmarkResult:
        t0 = time.time()
        ours = self.solve()
        elapsed = time.time() - t0

        ratios = None
        if self.baseline:
            ratios = {}
            for k in ("rms", "p2p_mean", "wire_length_m"):
                if k in ours and k in self.baseline:
                    a = float(ours[k]); b = float(self.baseline[k])
                    if b != 0:
                        ratios[k] = a / b

        result = BenchmarkResult(
            name=self.spec.name,
            reference=self.spec.reference,
            target_spec=asdict(self.spec),
            our_result=ours,
            published_baseline=self.baseline,
            ratio=ratios,
            elapsed_s=elapsed,
        )

        if json_out is not None:
            json_out = Path(json_out)
            json_out.write_text(
                json.dumps(asdict(result), indent=2, default=_json_default),
                encoding="utf-8",
            )
            print(f"  saved: {json_out}")
        return result


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer, np.floating)):
        return float(o)
    raise TypeError(f"Object of type {type(o)} is not JSON serializable")


def add_common_args(ap: argparse.ArgumentParser) -> None:
    """Add the standard CLI args every benchmark script uses."""
    ap.add_argument("--json", type=Path, default=None,
                    help="write the result as JSON to this file")
    ap.add_argument("--print-ratio", action="store_true",
                    help="print ratio table at the end if baseline known")
