#!/usr/bin/env python
"""Non-lattice hex ablation: does the HDiv-VIM Gram build stay practical on a hex
mesh that is NOT a translation lattice (so the translation-block cache cannot engage)?

The hex translation cache (m_hexUniformTransHosts) is auto-detected from geometry: it
turns ON only when every cell is a translated copy of a small template set (the
structured affine cube).  There is no env switch -- the ablation is DONE BY GEOMETRY.
This driver solves the SAME 1 m cube (uniform +z, linear mu_r=1000) at each N two ways:

  * LATTICE   : the affine structured hex map (translation cache ON)  -- Table 2's build
  * NONLATTICE: a smooth per-axis warp t -> t + a*sin(2*pi*t)/(2*pi) (a<1, so the map is
                monotone with positive Jacobian and fixes the [0,1] endpoints, keeping
                the cube domain) that makes every cell a DIFFERENT shape -> the cells are
                no longer translated copies -> the translation cache CANNOT engage, so the
                build runs on the geometry-agnostic techniques (i)-(iii) only.

Same N, same topology, same DOF count -- only the lattice-ness differs -- so the
LATTICE/NONLATTICE build-time ratio isolates the hex translation-cache speedup (the
companion of the wedge 8.9-14.9x already reported), and the NONLATTICE column is the
"(i)-(iii)-only on an arbitrary hex mesh" measurement the manuscript needs.

Benchmark Policy: timing -> mdx (idle-gated); one subprocess per case (per-case peak
working set).  Usage:
  python bench_hex_nonlattice.py --sizes 8,12,16,20          # driver (mdx)
  python bench_hex_nonlattice.py --case-n 12 --mode lattice  # one case, JSON on stdout
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import os
import platform
import subprocess
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
if len(_HERE.parents) > 3:
    _SRC = _HERE.parents[3] / "src"
    if _SRC.is_dir() and str(_SRC) not in sys.path:
        sys.path.insert(0, str(_SRC))

H0 = 200.0e3
MU_R = 1000.0
CUBE = 1.0
WARP_A = 0.35   # warp amplitude (<1 keeps the map monotone -> positive Jacobian)


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _peak_memory_mb() -> float:
    import psutil
    mem = psutil.Process(os.getpid()).memory_info()
    peak = getattr(mem, "peak_wset", None)
    return float(peak if peak is not None else mem.rss) / (1024.0 * 1024.0)


def _lattice_map(x, y, z):
    return (CUBE * (x - 0.5), CUBE * (y - 0.5), CUBE * (z - 0.5))


def _warp1(t):
    # t + a*sin(2 pi t)/(2 pi): fixes 0,1; derivative 1 + a*cos(2 pi t) > 0 for a < 1.
    return t + WARP_A * math.sin(2.0 * math.pi * t) / (2.0 * math.pi)


def _nonlattice_map(x, y, z):
    return (CUBE * (_warp1(x) - 0.5), CUBE * (_warp1(y) - 0.5), CUBE * (_warp1(z) - 0.5))


def run_case(n: int, mode: str) -> dict:
    """mode: 'lattice' (all caching on) / 'nonlattice' (warped: translation cache off by
    geometry, block ledger on) / 'naive' (lattice mesh + RADIA_HDIV_HEX_BLOCK_CACHE_LIMIT=1
    so both the translation block cache AND the block ledger are defeated -- every element-
    pair block is recomputed on each ACA entry request, the per-entry baseline)."""
    import numpy as np
    import ngsolve as ng
    from ngsolve.meshes import MakeStructured3DMesh
    import radia.vim as vim

    if mode == "naive":
        os.environ["RADIA_HDIV_HEX_BLOCK_CACHE_LIMIT"] = "1"
    mapping = _nonlattice_map if mode == "nonlattice" else _lattice_map
    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n, mapping=mapping)
        res = vim.Solve(mesh, MU_R, ng.CoefficientFunction((0.0, 0.0, H0)), order=1)

    return dict(
        N=int(n), mode=mode, ndof=int(res["ndof"]), n_el=int(res["n_el"]),
        n_charge=int(res["n_charge"]),
        charge_gram_wall_s=float(res["charge_gram_wall_s"]),
        iters=int(res["iters"]),
        M_avg_z=float(np.asarray(res["M_avg"], float)[2]),
        peak_memory_mb=_peak_memory_mb(),
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sizes", default="8,12,16,20",
                    help="N for the lattice-vs-nonlattice (translation-cache) sweep")
    ap.add_argument("--naive-sizes", default="6,8,10",
                    help="N for the naive baseline (block ledger off ~80x, so small N only)")
    ap.add_argument("--case-n", type=int, default=None)
    ap.add_argument("--mode", choices=["lattice", "nonlattice", "naive"], default="lattice")
    ap.add_argument("--out", default=str(_HERE.parent / "results_hex_nonlattice.json"))
    args = ap.parse_args()

    if args.case_n is not None:
        print(json.dumps(run_case(args.case_n, args.mode)))
        return 0

    def _case(n, mode):
        proc = subprocess.run([sys.executable, str(_HERE), "--case-n", str(n), "--mode", mode],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"case N={n} {mode} failed:\n{proc.stdout}\n{proc.stderr}")
        return json.loads(proc.stdout.strip().splitlines()[-1])

    # Part 1: translation-cache ablation (lattice vs non-lattice) + the non-lattice hex
    # build numbers, at the Table-2 sizes.
    trans = []
    for n in [int(s) for s in args.sizes.split(",") if s.strip()]:
        lat, non = _case(n, "lattice"), _case(n, "nonlattice")
        ratio = non["charge_gram_wall_s"] / lat["charge_gram_wall_s"]
        reldiff = abs(non["M_avg_z"] - lat["M_avg_z"]) / abs(lat["M_avg_z"])
        trans.append(dict(N=n, ndof=lat["ndof"], lattice=lat, nonlattice=non,
                          nonlattice_over_lattice=ratio, M_avg_z_rel_diff=reldiff))
        print(f"[trans] N={n:2d} ndof={lat['ndof']:7d}  lattice={lat['charge_gram_wall_s']:8.2f}s "
              f"nonlattice={non['charge_gram_wall_s']:8.2f}s  (x{ratio:.2f})  "
              f"M_avg reldiff={reldiff:.2e}", flush=True)

    # Part 2: naive baseline (all block/translation caching defeated) vs full, small N.
    naive = []
    for n in [int(s) for s in args.naive_sizes.split(",") if s.strip()]:
        full, nv = _case(n, "lattice"), _case(n, "naive")
        ratio = nv["charge_gram_wall_s"] / full["charge_gram_wall_s"]
        naive.append(dict(N=n, ndof=full["ndof"], full=full, naive=nv, naive_over_full=ratio))
        print(f"[naive] N={n:2d} ndof={full['ndof']:7d}  full={full['charge_gram_wall_s']:8.2f}s "
              f"naive={nv['charge_gram_wall_s']:9.2f}s  (x{ratio:.1f})", flush=True)

    data = dict(
        timestamp=_now(), hostname=platform.node(),
        benchmark="hex_gram_build_ablation",
        problem=dict(cube_size_m=CUBE, H0_A_per_m=H0, mu_r=MU_R,
                     warp_amplitude=WARP_A, mesh="structured hex NxNxN", order="RT1"),
        translation_cache=trans,
        naive_baseline=naive,
    )
    Path(args.out).write_text(json.dumps(data, indent=1))
    print(f"[hex_nonlattice] wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
