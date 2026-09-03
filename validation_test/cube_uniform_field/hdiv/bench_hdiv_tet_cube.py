#!/usr/bin/env python
"""Tet-mesh HDiv-VIM cube benchmark for the soft-iron uniform-field problem.

This is the unstructured-tet counterpart of bench_hdiv_cube.py.  It uses the
current production path, radia.vim.Solve, not the retired Radia object/MMPM
benchmark. Useful sizes are validation-class workloads and should run on
hibino first, or on mdx only when hibino is unavailable and its CI queue is
idle.
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


_HERE = Path(__file__).resolve()
REPO_ROOT = _HERE.parents[3] if len(_HERE.parents) > 3 else _HERE.parent
_SRC = REPO_ROOT / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

H_EXT = 200_000.0
CUBE_SIZE = 1.0
DEFAULT_MU_R = 1000.0

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


def _make_tet_cube(maxh: float, size: float):
    import ngsolve as ng
    from netgen.csg import CSGeometry, OrthoBrick, Pnt

    half = 0.5 * float(size)
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-half, -half, -half), Pnt(half, half, half)))
    return ng.Mesh(geo.GenerateMesh(maxh=float(maxh)))


def _solve(mesh, solve_kwargs: dict[str, Any]):
    import ngsolve as ng
    from radia.vim import Solve

    with ng.TaskManager():
        return Solve(mesh, **solve_kwargs)


def _case_name(material: str, maxh: float) -> str:
    h = ("%0.6f" % float(maxh)).rstrip("0").rstrip(".").replace(".", "p")
    return f"{material}_tet_maxh{h}"


def run_case(args: argparse.Namespace, maxh: float, material: str, bh_table: list[list[float]]) -> dict[str, Any]:
    import ngsolve as ng
    import radia

    mesh = _make_tet_cube(maxh, args.cube_size)
    solve_kwargs: dict[str, Any] = {
        "H_ext": ng.CoefficientFunction((0.0, 0.0, args.h_ext)),
        "order": args.order,
        "gram_eps": args.gram_eps,
        "leaf": args.leaf,
        "eta": args.eta,
        "tol": args.tol,
        "maxit": args.maxit,
        "nl_tol": args.nl_tol,
        "nl_maxit": args.nl_maxit,
        "linear_solver": args.linear_solver,
        "preconditioner": args.preconditioner,
        "newton_inner_tol": args.newton_inner_tol,
        "newton_warmstart": args.newton_warmstart,
        "newton_continuation": args.newton_continuation,
        "newton_reuse_tangent_steps": args.newton_reuse_tangent_steps,
        "newton_cg_x0": args.newton_cg_x0,
    }
    if material == "linear":
        solve_kwargs["mu_r"] = args.mu_r
    elif material == "nonlinear":
        solve_kwargs["bh_table"] = bh_table
        solve_kwargs["nonlinear_solver"] = args.nonlinear_solver
    else:
        raise ValueError(f"unknown material: {material}")

    t0 = time.perf_counter()
    res = _solve(mesh, solve_kwargs)
    wall_s = time.perf_counter() - t0

    m_avg = [float(x) for x in np.asarray(res["M_avg"], dtype=float)]
    mz = m_avg[2]
    transverse_ratio = math.hypot(m_avg[0], m_avg[1]) / abs(mz) if mz else None
    row: dict[str, Any] = {
        "case": _case_name(material, maxh),
        "material": material,
        "mesh_kind": "tet",
        "maxh": float(maxh),
        "nonlinear_profile": args.nonlinear_profile,
        "newton_cg_x0": bool(args.newton_cg_x0),
        "n_el": int(res["n_el"]),
        "ndof": int(res["ndof"]),
        "n_charge": int(res["n_charge"]),
        "iters": int(res["iters"]),
        "demag": float(res["demag"]),
        "M_avg_A_per_m": m_avg,
        "M_avg_z_A_per_m": mz,
        "transverse_ratio": transverse_ratio,
        "wall_s": wall_s,
        "setup_wall_s": res.get("setup_wall_s"),
        "fes_wall_s": res.get("fes_wall_s"),
        "charge_gram_wall_s": res.get("charge_gram_wall_s"),
        "charge_basis_wall_s": res.get("charge_basis_wall_s"),
        "charge_gram_cpp_wall_s": res.get("charge_gram_cpp_wall_s"),
        "projection_wall_s": res.get("projection_wall_s"),
        "demag_probe_wall_s": res.get("demag_probe_wall_s"),
        "solve_wall_s": res.get("solve_wall_s"),
        "cpp_solve_timings": res.get("cpp_solve_timings"),
        "nonlinear_solve_stats": res.get("nonlinear_solve_stats"),
        "post_wall_s": res.get("post_wall_s"),
        "total_wall_s_internal": res.get("total_wall_s_internal"),
        "linear_solver": res.get("linear_solver"),
        "preconditioner": res.get("preconditioner"),
        "preconditioner_requested": res.get("preconditioner_requested"),
        "preconditioner_policy": res.get("preconditioner_policy"),
        "gram_backend": res.get("gram_backend"),
        "order": int(res.get("order", args.order)),
        "nonlinear": bool(res.get("nonlinear", False)),
        "hmat_stats": res.get("hmat_stats"),
        "memory": _memory_mb(),
        "mesh_materials": list(mesh.GetMaterials()),
        "mesh_boundaries": list(mesh.GetBoundaries()),
        "radia_version": getattr(radia, "__version__", None),
    }
    for key in [
        "solve_total_s",
        "solve_factor_s",
        "solve_prec_s",
        "solve_bx_s",
        "solve_gmatvec_s",
        "solve_btx_s",
        "solve_mass_s",
        "solve_dot_s",
        "solve_ax_total_s",
        "solve_ax_other_s",
        "solve_pcg_update_s",
        "solve_apply_count",
        "solve_prec_count",
        "solve_dot_count",
        "hmatvec_total_s",
        "hmatvec_zero_s",
        "hmatvec_permute_s",
        "hmatvec_leaf_s",
        "hmatvec_reduce_s",
        "hmatvec_meta_s",
        "hmatvec_lowrank_flop_est",
        "hmatvec_dense_flop_est",
        "hmatvec_calls",
        "hmatvec_lowrank_leaves",
        "hmatvec_dense_leaves",
        "hmatvec_mirrored_upper_leaves",
        "hmatvec_diagonal_leaves",
        "hmatvec_skipped_lower_leaves",
        "hmatvec_last_nd",
        "hmatvec_last_nthr",
        "nonlinear_newton_iters",
        "nonlinear_warmstart_solves",
        "nonlinear_linear_inner_iters",
        "nonlinear_line_search_backtracks",
        "nonlinear_tangent_assemblies",
        "nonlinear_tangent_reuses",
        "nonlinear_fresh_tangent_retries",
        "nonlinear_final_rel_step",
        "nonlinear_final_settled_iters",
    ]:
        if key in res:
            row[key] = res.get(key)
    if material == "linear":
        chi = args.mu_r - 1.0
        ref = chi * args.h_ext / (1.0 + chi / 3.0)
        row["uniform_demag_approx_Mz_A_per_m"] = ref
        row["uniform_demag_approx_rel_difference"] = abs(mz - ref) / abs(ref)
    return row


def _materials(selection: str) -> list[str]:
    if selection == "both":
        return ["linear", "nonlinear"]
    return [selection]


def _apply_nonlinear_profile(args: argparse.Namespace) -> None:
    """Resolve nonlinear timing defaults after argparse.

    `strict` keeps solver-oriented defaults.  `mdx-scaling` is the timing lane
    used for the 2026-07-09 nonlinear preconditioner sweep; explicit CLI values
    still win because these fields default to None in argparse.
    """
    if args.preconditioner is None:
        args.preconditioner = "auto"
    if args.nl_tol is None:
        args.nl_tol = 3.0e-4 if args.nonlinear_profile == "mdx-scaling" else 1.0e-6
    if args.nl_maxit is None:
        args.nl_maxit = 300
    if args.newton_inner_tol is None:
        args.newton_inner_tol = "auto"
    if args.newton_warmstart is None:
        args.newton_warmstart = "linear"
    if args.newton_continuation is None:
        args.newton_continuation = 2 if args.nonlinear_profile == "mdx-scaling" else 1
    if args.newton_reuse_tangent_steps is None:
        args.newton_reuse_tangent_steps = 3 if args.nonlinear_profile == "mdx-scaling" else 1


def _build_payload(args: argparse.Namespace, rows: list[dict[str, Any]], bh_source: str) -> dict[str, Any]:
    return {
        "generated_at_utc": _now_utc(),
        "driver": str(Path(__file__).resolve()),
        "source_manuscript": args.source_manuscript,
        "presentation_dir": args.presentation_dir,
        "problem": {
            "geometry": "1 m cube centered at origin",
            "mesh": "Netgen CSG unstructured tet cube, maxh-controlled",
            "method": "BDM1 HDiv-VIM via radia.vim.Solve",
            "H_ext_A_per_m": args.h_ext,
            "mu_r_linear": args.mu_r,
            "bh_table_source": bh_source,
        },
        "parameters": {
            "maxh_values": args.maxh_values,
            "material": args.material,
            "order": args.order,
            "gram_eps": args.gram_eps,
            "leaf": args.leaf,
            "eta": args.eta,
            "tol": args.tol,
            "maxit": args.maxit,
            "nl_tol": args.nl_tol,
            "nl_maxit": args.nl_maxit,
            "linear_solver": args.linear_solver,
            "preconditioner": args.preconditioner,
            "nonlinear_solver": args.nonlinear_solver,
            "nonlinear_profile": args.nonlinear_profile,
            "newton_inner_tol": args.newton_inner_tol,
            "newton_warmstart": args.newton_warmstart,
            "newton_continuation": args.newton_continuation,
            "newton_reuse_tangent_steps": args.newton_reuse_tangent_steps,
            "newton_cg_x0": args.newton_cg_x0,
        },
        "runtime": {
            "python": sys.version,
            "platform": platform.platform(),
            "machine": platform.node(),
        },
        "results": rows,
    }


def _write_payload(path: str, payload: dict[str, Any]) -> None:
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--maxh-values", nargs="+", type=float, default=[0.4, 0.3])
    ap.add_argument("--material", choices=["linear", "nonlinear", "both"], default="both")
    ap.add_argument(
        "--output",
        type=str,
        default=str(_HERE.parent / "results_hdiv_tet_cube.json"),
        help="JSON result path (default: validation directory)",
    )
    ap.add_argument("--bh-table", type=str, default="")
    ap.add_argument("--source-manuscript", type=str, default="")
    ap.add_argument("--presentation-dir", type=str, default="")
    ap.add_argument("--cube-size", type=float, default=CUBE_SIZE)
    ap.add_argument("--h-ext", type=float, default=H_EXT)
    ap.add_argument("--mu-r", type=float, default=DEFAULT_MU_R)
    ap.add_argument("--order", type=int, choices=[1, 2], default=1)
    ap.add_argument("--gram-eps", type=float, default=1.0e-4)
    ap.add_argument("--leaf", type=int, default=32)
    ap.add_argument("--eta", type=float, default=2.0)
    ap.add_argument("--tol", type=float, default=1.0e-8)
    ap.add_argument("--maxit", type=int, default=4000)
    ap.add_argument("--nl-tol", type=float, default=None)
    ap.add_argument("--nl-maxit", type=int, default=None)
    ap.add_argument("--linear-solver", choices=["auto", "cpp-cg", "gmres"], default="auto")
    ap.add_argument("--preconditioner", choices=["auto", "mass-riesz", "jacobi"], default=None)
    ap.add_argument(
        "--nonlinear-solver",
        choices=["energy-newton", "picard-mass-riesz", "picard-energy"],
        default="energy-newton",
    )
    ap.add_argument("--nonlinear-profile", choices=["strict", "mdx-scaling"], default="strict",
                    help="Resolve nonlinear defaults. mdx-scaling uses nl_tol=3e-4, continuation=2, "
                         "and tangent reuse=3 unless explicitly overridden.")
    ap.add_argument("--newton-inner-tol", default=None,
                    help="'auto', 'fixed', or a numeric inner CG tolerance floor for energy-newton.")
    ap.add_argument("--newton-warmstart", choices=["linear", "picard", "none"], default=None)
    ap.add_argument("--newton-continuation", type=int, default=None)
    ap.add_argument("--newton-reuse-tangent-steps", type=int, default=None)
    ap.add_argument("--newton-cg-x0", action="store_true",
                    help="Use previous Newton direction as the inner CG initial guess.")
    args = ap.parse_args()
    _apply_nonlinear_profile(args)

    bh_table, bh_source = _load_bh_table(args.bh_table or None)
    rows = []
    for maxh in args.maxh_values:
        for material in _materials(args.material):
            print(f"[hdiv-tet-cube] {material} maxh={maxh}", flush=True)
            rows.append(run_case(args, maxh, material, bh_table))
            if args.output:
                _write_payload(args.output, _build_payload(args, rows, bh_source))

    payload = _build_payload(args, rows, bh_source)
    text = json.dumps(payload, ensure_ascii=False, indent=2)
    if args.output:
        _write_payload(args.output, payload)
        print(f"[hdiv-tet-cube] wrote {Path(args.output)}")
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
