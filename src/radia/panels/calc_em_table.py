"""
Offline builder for the 2D EM table (|H_t|, T) -> Z_s.

Replaces the outer Karl FEM iteration in calc_fem_kelvin.py
(``--impedance esim``) for the engineering use case where the EM
problem is quasi-static at the heat timescale.  Pre-computes the
effective surface impedance Z_s on a 2D grid spanned by the tangential
field magnitude H_t [A/m] and the workpiece temperature T [degC].  The
runtime heat solver (``calc_heat_with_em_table.py``) bilinearly
interpolates this table per surface DOF, per timestep, instead of
re-solving the EM problem.

CLI
---
    python calc_em_table.py \
        --material steel --frequency 50e3 \
        --H-grid-min 1e2  --H-grid-max 1e5 --n-H 20 \
        --T-grid-min 20   --T-grid-max 800 --n-T 20 \
        --sigma-T-curve sigma_T_steel.csv \
        --output em_table_steel_50kHz.npz

``--sigma-T-curve`` is a 2-column CSV (``T_celsius, sigma_S_per_m``)
that overrides the constant ``EMMaterial.sigma``.  Outside the
tabulated T range the curve is clamped (no extrapolation).  When
omitted, the material's preset constant sigma is used at all
temperatures (table is then "1D in H, replicated across T") -- useful
for sanity checks against the existing 1D Karl path at fixed T.

BH curve is treated as T-invariant in v1; T-dependent BH is on the
roadmap (see CLAUDE.md "Production roadmap" note).  This is the
single dominant simplification of the v1 table; the dominant
T-dependence in IH steel below the Curie point is sigma, not BH.

Output ``.npz`` schema
----------------------
    H_grid    (n_H,)        float, A/m, log-spaced
    T_grid    (n_T,)        float, degC, linear
    Zs_re     (n_H, n_T)    float, Ohm
    Zs_im     (n_H, n_T)    float, Ohm
    q_surf    (n_H, n_T)    float, W/m^2   (= 0.5 * Re(Z_s) * H^2)
    meta      0-d str       JSON: material, frequency, sigma_T_path,
                                  n_nodes, num_skin_depths, build_time_s
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

_panels_dir = os.path.dirname(os.path.abspath(__file__))
if _panels_dir not in sys.path:
    sys.path.insert(0, _panels_dir)

from calc_common import setup_paths, progress, calc_main  # noqa: E402

setup_paths()
from em_material import EMMaterial  # noqa: E402
from esim_cell_problem import ESIMCellProblemSolver  # noqa: E402


def _log(msg):
    progress("EM_TABLE", msg)


def load_sigma_T_csv(path):
    """Load a 2-column ``T_celsius, sigma_S_per_m`` CSV.

    Returns a callable ``sigma_at(T_celsius)`` that linearly
    interpolates and CLAMPS outside the tabulated range (no
    extrapolation -- pinning sigma at the endpoints is the only
    physically defensible behaviour past the data).
    """
    data = np.loadtxt(path, delimiter=",", comments="#")
    if data.ndim != 2 or data.shape[1] < 2:
        raise ValueError(
            f"--sigma-T-curve must be 2-col CSV (T_celsius, sigma_S_per_m), "
            f"got shape {data.shape} from {path}")
    T_arr = data[:, 0].astype(float)
    s_arr = data[:, 1].astype(float)
    order = np.argsort(T_arr)
    T_arr = T_arr[order]
    s_arr = s_arr[order]
    if np.any(s_arr <= 0):
        raise ValueError(
            f"--sigma-T-curve contains non-positive sigma values: {s_arr}")

    def sigma_at(T_c):
        return float(np.interp(T_c, T_arr, s_arr,
                               left=s_arr[0], right=s_arr[-1]))

    sigma_at.T_min = float(T_arr[0])
    sigma_at.T_max = float(T_arr[-1])
    return sigma_at


def build_table(material, frequency, H_grid, T_grid, sigma_at,
                n_nodes=100, num_skin_depths=10):
    """Sweep (H, T) and return (Zs_re, Zs_im, q_surf) arrays.

    Each cell instantiates a fresh ESIMCellProblemSolver (cheap; the
    mesh build is dominated by a few hundred floats) so the sigma at
    that T can flow into the rho = 1/sigma the solver bakes into its
    PDE coefficients.
    """
    n_H = len(H_grid)
    n_T = len(T_grid)
    Zs_re = np.zeros((n_H, n_T), dtype=float)
    Zs_im = np.zeros((n_H, n_T), dtype=float)

    bh = material.bh_curve
    t0 = time.perf_counter()
    for j, T_c in enumerate(T_grid):
        sigma_T = sigma_at(T_c)
        solver = ESIMCellProblemSolver(
            bh_curve=bh, sigma=sigma_T, frequency=frequency,
            n_nodes=n_nodes, num_skin_depths=num_skin_depths)
        for i, H in enumerate(H_grid):
            res = solver.solve(float(H))
            Z = complex(res["Z"])
            Zs_re[i, j] = Z.real
            Zs_im[i, j] = Z.imag
        _log(f"col {j+1}/{n_T} T={T_c:.1f}C sigma={sigma_T:.3e} "
             f"done t={time.perf_counter() - t0:.1f}s")

    q_surf = 0.5 * Zs_re * (H_grid[:, None] ** 2)
    return Zs_re, Zs_im, q_surf


def run(args):
    setup_paths()
    if args.H_grid_min <= 0:
        raise ValueError("--H-grid-min must be positive (log spacing)")
    if args.H_grid_max <= args.H_grid_min:
        raise ValueError("--H-grid-max must exceed --H-grid-min")
    if args.T_grid_max < args.T_grid_min:
        raise ValueError("--T-grid-max must exceed --T-grid-min")
    if args.n_H < 2 or args.n_T < 2:
        raise ValueError("--n-H and --n-T must each be >= 2")

    material = EMMaterial.from_args(args)
    _log(f"material={material.name} sigma_default={material.sigma:.3e} "
         f"mu_r={material.mu_r:g} bh={'yes' if material.bh_curve else 'no'}")

    if args.sigma_T_curve:
        sigma_at = load_sigma_T_csv(args.sigma_T_curve)
        _log(f"sigma(T) CSV: {args.sigma_T_curve} "
             f"range T=[{sigma_at.T_min:.1f},{sigma_at.T_max:.1f}]C")
        sigma_T_path = os.path.abspath(args.sigma_T_curve)
    else:
        sigma_const = material.sigma
        sigma_at = lambda T_c: sigma_const  # noqa: E731
        sigma_at.T_min = float("-inf")
        sigma_at.T_max = float("+inf")
        _log(f"sigma(T) = constant {sigma_const:.3e} "
             f"(no --sigma-T-curve given; table T axis is redundant)")
        sigma_T_path = ""

    H_grid = np.logspace(np.log10(args.H_grid_min),
                         np.log10(args.H_grid_max), args.n_H)
    T_grid = np.linspace(args.T_grid_min, args.T_grid_max, args.n_T)
    _log(f"H_grid: {args.n_H} pts [{H_grid[0]:.2e}, {H_grid[-1]:.2e}] A/m")
    _log(f"T_grid: {args.n_T} pts [{T_grid[0]:.1f}, {T_grid[-1]:.1f}] C")

    t0 = time.perf_counter()
    Zs_re, Zs_im, q_surf = build_table(
        material, args.frequency, H_grid, T_grid, sigma_at,
        n_nodes=args.n_nodes, num_skin_depths=args.num_skin_depths)
    build_time_s = time.perf_counter() - t0

    out_path = args.em_table_output
    if not out_path:
        raise ValueError("--em-table-output is required")
    out_path = os.path.abspath(out_path)
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    meta = {
        "material": material.name,
        "frequency": args.frequency,
        "sigma_T_curve": sigma_T_path,
        "sigma_default": material.sigma,
        "mu_r": material.mu_r,
        "n_nodes": args.n_nodes,
        "num_skin_depths": args.num_skin_depths,
        "build_time_s": build_time_s,
        "schema_version": 1,
    }
    np.savez(out_path,
             H_grid=H_grid, T_grid=T_grid,
             Zs_re=Zs_re, Zs_im=Zs_im, q_surf=q_surf,
             meta=np.asarray(json.dumps(meta)))
    _log(f"saved {out_path} ({os.path.getsize(out_path)/1024:.1f} KB)")

    return {
        "em_table": out_path,
        "n_H": int(args.n_H),
        "n_T": int(args.n_T),
        "frequency": args.frequency,
        "material": material.name,
        "Zs_re_min": float(Zs_re.min()),
        "Zs_re_max": float(Zs_re.max()),
        "q_surf_min": float(q_surf.min()),
        "q_surf_max": float(q_surf.max()),
        "build_time_s": build_time_s,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Build (|H_t|, T) -> Z_s table for the engineering "
                    "heat-only IH workflow.")
    # Material (shared --material / --sigma / --mu-r / --bh-file).
    from em_material import add_material_args
    add_material_args(parser, default_material="steel")

    parser.add_argument("--frequency", type=float, required=True,
                        help="Operating frequency [Hz]")
    parser.add_argument("--H-grid-min", type=float, default=1e2,
                        help="Min |H_t| [A/m] (must be > 0)")
    parser.add_argument("--H-grid-max", type=float, default=1e5,
                        help="Max |H_t| [A/m]")
    parser.add_argument("--n-H", type=int, default=20,
                        help="Number of H points (log-spaced)")
    parser.add_argument("--T-grid-min", type=float, default=20.0,
                        help="Min temperature [degC]")
    parser.add_argument("--T-grid-max", type=float, default=800.0,
                        help="Max temperature [degC]")
    parser.add_argument("--n-T", type=int, default=20,
                        help="Number of T points (linear)")
    parser.add_argument("--sigma-T-curve", default="",
                        help="2-col CSV (T_celsius, sigma_S_per_m).  "
                             "When given, overrides EMMaterial.sigma at "
                             "each grid T (linear interp, clamped).")
    parser.add_argument("--n-nodes", type=int, default=100,
                        help="ESIM cell-problem 1D mesh nodes")
    parser.add_argument("--num-skin-depths", type=float, default=10.0,
                        help="ESIM domain depth in skin depths")
    parser.add_argument("--em-table-output", default="",
                        help=".npz output path (required)")
    calc_main(run, parser)


if __name__ == "__main__":
    main()
