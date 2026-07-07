#!/usr/bin/env python
"""HDiv-VIM cube benchmark driver for the soft-iron uniform-field problem.

This is the HDiv counterpart of the MMPM cube benchmark used in the
2026-08-25 static/rotating-machine manuscript:

  * 1 m x 1 m x 1 m cube, centered at the origin
  * structured pure-hex mesh, N x N x N
  * applied field H0 = 200 kA/m in +z
  * RT1 HDiv-VIM solve through radia.vim.Solve
  * HACApK charge-Gram stats recorded from the solver result

The script is intentionally in validation_test, not tests: useful sizes are
solver-heavy and should be run on mdx for publication-grade timing.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))

MU0 = 4.0e-7 * math.pi
H_EXT = 200_000.0
CUBE_SIZE = 1.0
DEFAULT_MU_R = 1000.0

# Small built-in fallback.  For publication runs, pass --bh-table pointing at
# the BH.txt used by the manuscript so the material law is identical.
FALLBACK_BH_TABLE = [
    [0.0, 0.0],
    [100.0, 0.1],
    [200.0, 0.3],
    [500.0, 0.8],
    [1000.0, 1.2],
    [2000.0, 1.5],
    [5000.0, 1.7],
    [10000.0, 1.8],
    [50000.0, 2.0],
    [100000.0, 2.1],
    [318000.0, 2.61],
]


def _now_utc() -> str:
    return _dt.datetime.now(_dt.UTC).replace(microsecond=0).isoformat()


def _load_bh_table(path: str | None) -> tuple[list[list[float]], str]:
    if not path:
        return FALLBACK_BH_TABLE, "built-in fallback BH table"
    p = Path(path)
    rows: list[list[float]] = []
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        if len(parts) < 2:
            continue
        rows.append([float(parts[0]), float(parts[1])])
    if len(rows) < 3:
        raise ValueError(f"BH table too short: {p}")
    return rows, str(p)


def _memory_mb() -> dict[str, float | None]:
    try:
        import psutil
    except ImportError:
        return {"rss_mb": None, "peak_wset_mb": None}
    proc = psutil.Process(os.getpid())
    info = proc.memory_info()
    rss = info.rss / (1024.0 * 1024.0)
    peak = getattr(info, "peak_wset", None)
    return {
        "rss_mb": rss,
        "peak_wset_mb": None if peak is None else peak / (1024.0 * 1024.0),
    }


def _make_hex_cube(n: int, size: float):
    from ngsolve.meshes import MakeStructured3DMesh

    return MakeStructured3DMesh(
        hexes=True,
        nx=n,
        ny=n,
        nz=n,
        mapping=lambda x, y, z: (
            size * (x - 0.5),
            size * (y - 0.5),
            size * (z - 0.5),
        ),
    )


def _solve_with_retry(mesh, solve_kwargs: dict[str, Any], retries: int):
    import ngsolve as ng
    from radia.vim import Solve

    last: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with ng.TaskManager():
                return Solve(mesh, **solve_kwargs), attempt
        except RuntimeError as exc:
            if "GetTrafo lattice evaluation unstable" not in str(exc):
                raise
            last = exc
    raise RuntimeError(f"vim.Solve failed after {retries} retries: {last}") from last


def run_case(args: argparse.Namespace, n: int, material: str, bh_table: list[list[float]]) -> dict[str, Any]:
    import ngsolve as ng
    import radia

    mesh = _make_hex_cube(n, args.cube_size)
    solve_kwargs: dict[str, Any] = {
        "H_ext": ng.CoefficientFunction((0.0, 0.0, args.h_ext)),
        "gram_eps": args.gram_eps,
        "leaf": args.leaf,
        "eta": args.eta,
        "tol": args.tol,
        "maxit": args.maxit,
        "nl_tol": args.nl_tol,
        "nl_maxit": args.nl_maxit,
    }
    if material == "linear":
        solve_kwargs["mu_r"] = args.mu_r
    elif material == "nonlinear":
        solve_kwargs["bh_table"] = bh_table
    else:
        raise ValueError(f"unknown material: {material}")

    t0 = time.perf_counter()
    res, attempts = _solve_with_retry(mesh, solve_kwargs, args.retries)
    wall_s = time.perf_counter() - t0

    m_avg = [float(x) for x in np.asarray(res["M_avg"], dtype=float)]
    mz = m_avg[2]
    transverse_ratio = (
        math.hypot(m_avg[0], m_avg[1]) / abs(mz)
        if mz else None
    )
    row: dict[str, Any] = {
        "case": f"{material}_N{n}",
        "material": material,
        "N": int(n),
        "n_el": int(res["n_el"]),
        "ndof": int(res["ndof"]),
        "n_charge": int(res["n_charge"]),
        "iters": int(res["iters"]),
        "demag": float(res["demag"]),
        "M_avg_A_per_m": m_avg,
        "M_avg_z_A_per_m": mz,
        "transverse_ratio": transverse_ratio,
        "wall_s": wall_s,
        "attempts": attempts,
        "linear_solver": res.get("linear_solver"),
        "gram_backend": res.get("gram_backend"),
        "order": int(res.get("order", 1)),
        "nonlinear": bool(res.get("nonlinear", False)),
        "hmat_stats": res.get("hmat_stats"),
        "memory": _memory_mb(),
        "radia_version": getattr(radia, "__version__", None),
    }
    if material == "linear":
        chi = args.mu_r - 1.0
        # This is a useful low-susceptibility sanity estimate, not a high-mu
        # cube truth value: the HDiv solution permits nonuniform M, whereas the
        # closed expression assumes a uniform-magnetization demag factor.
        ref = chi * args.h_ext / (1.0 + chi / 3.0)
        row["uniform_demag_approx_Mz_A_per_m"] = ref
        row["uniform_demag_approx_rel_difference"] = abs(mz - ref) / abs(ref)
    return row


def _materials(selection: str) -> list[str]:
    if selection == "both":
        return ["linear", "nonlinear"]
    return [selection]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sizes", nargs="+", type=int, default=[2, 4])
    ap.add_argument("--material", choices=["linear", "nonlinear", "both"], default="both")
    ap.add_argument("--output", type=str, default="")
    ap.add_argument("--bh-table", type=str, default="")
    ap.add_argument("--source-manuscript", type=str, default="")
    ap.add_argument("--presentation-dir", type=str, default="")
    ap.add_argument("--cube-size", type=float, default=CUBE_SIZE)
    ap.add_argument("--h-ext", type=float, default=H_EXT)
    ap.add_argument("--mu-r", type=float, default=DEFAULT_MU_R)
    ap.add_argument("--gram-eps", type=float, default=1.0e-4)
    ap.add_argument("--leaf", type=int, default=32)
    ap.add_argument("--eta", type=float, default=2.0)
    ap.add_argument("--tol", type=float, default=1.0e-8)
    ap.add_argument("--maxit", type=int, default=4000)
    ap.add_argument("--nl-tol", type=float, default=1.0e-6)
    ap.add_argument("--nl-maxit", type=int, default=300)
    ap.add_argument("--retries", type=int, default=5)
    args = ap.parse_args()

    bh_table, bh_source = _load_bh_table(args.bh_table or None)
    rows = []
    for n in args.sizes:
        for material in _materials(args.material):
            print(f"[hdiv-cube] {material} N={n}", flush=True)
            rows.append(run_case(args, n, material, bh_table))

    payload = {
        "generated_at_utc": _now_utc(),
        "driver": str(Path(__file__).resolve()),
        "source_manuscript": args.source_manuscript,
        "presentation_dir": args.presentation_dir,
        "problem": {
            "geometry": "1 m cube centered at origin",
            "mesh": "structured pure hex N x N x N",
            "method": "RT1 HDiv-VIM via radia.vim.Solve",
            "H_ext_A_per_m": args.h_ext,
            "mu_r_linear": args.mu_r,
            "bh_table_source": bh_source,
        },
        "parameters": {
            "sizes": args.sizes,
            "material": args.material,
            "gram_eps": args.gram_eps,
            "leaf": args.leaf,
            "eta": args.eta,
            "tol": args.tol,
            "maxit": args.maxit,
            "nl_tol": args.nl_tol,
            "nl_maxit": args.nl_maxit,
        },
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.node(),
        },
        "results": rows,
    }
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if args.output:
        out = Path(args.output)
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(text + "\n", encoding="utf-8")
        print(f"[hdiv-cube] wrote {out}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
