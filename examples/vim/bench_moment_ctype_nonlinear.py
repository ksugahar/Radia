"""Nonlinear C-yoke gap-field benchmark for the collocation-MMMM coarse tier (method 2).

The engineering question: the GAP-CENTER field of a voxelized C-yoke with a NONLINEAR
BH curve (MatSatIsoTab), at ~165k DoF (nside=54 -> 27,720 hex -> 166,320 DoF), solved by
the production HACApK-BiCGSTAB moment route (analytic kernel + Anderson(1) Picard, both
default since radia 4.95.5).

Geometry = the same voxelized C-yoke as bench_moment_solvers.py (outer 0.12 m square,
0.07 m cavity, gap strip removed for x>0.018 so the two pole arms face across y):
uniform in-plane y-drive B_ext = mu0*H0 excites the C loop; gap center at [0.045,0,0].

Self-calibrating drive: a small nside (default 24, 12.7k DoF) is solved at ascending H0
until max element |M| crosses M_KNEE (=9e5 A/m ~ B_local 1.2 T, the BH knee), so the big
case runs in the genuinely NONLINEAR partial-saturation regime.  A LINEAR reference
(mu_r = the BH initial slope) runs at the same H0 and size to quantify the saturation
effect on the gap field.

Self-contained (no repo needed on the bench machine).  JSON per Benchmark Policy.
"""
import argparse
import json
import os
import platform
import time
from datetime import datetime

import numpy as np
import psutil
import radia as rad

MU0 = 4.0e-7 * np.pi

# Dense monotone steel-like curve (knee ~1.2-1.5 T, Bs ~2 T); initial slope
# mu_r ~ 0.1/(mu0*100) ~ 796.  MatSatIsoTab interpolates [[H, B], ...].
BH_DATA = [[0.0, 0.0], [50.0, 0.05], [100.0, 0.1], [200.0, 0.25], [400.0, 0.6],
           [600.0, 0.9], [800.0, 1.05], [1000.0, 1.2], [2000.0, 1.45],
           [5000.0, 1.65], [10000.0, 1.75], [20000.0, 1.85], [50000.0, 2.0]]
MU_R_LINEAR_REF = 0.1 / (MU0 * 100.0)      # BH initial slope as the linear reference
M_KNEE = 9.0e5                              # A/m; max|M| above this = partial saturation reached

GAP_PROBE = [0.045, 0.0, 0.0]
EXT_PROBES = [[0.09, 0.03, 0.005], [0.0, 0.09, 0.01]]


def _inside_cyoke(cx, cy):
    if not (-0.06 <= cx <= 0.06 and -0.06 <= cy <= 0.06):
        return False
    if -0.035 <= cx <= 0.035 and -0.035 <= cy <= 0.035:
        return False
    if cx > 0.018 and -0.035 <= cy <= 0.035:
        return False
    return True


def build_ctype(nside, mat_factory):
    nxy, nz = nside, max(2, nside // 3)
    xs = np.linspace(-0.06, 0.06, nxy + 1)
    zs = np.linspace(-0.02, 0.02, nz + 1)
    objs = []
    for k in range(nz):
        for j in range(nxy):
            for i in range(nxy):
                if not _inside_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])):
                    continue
                x0, x1, y0, y1 = xs[i], xs[i + 1], xs[j], xs[j + 1]
                z0, z1 = zs[k], zs[k + 1]
                v = [[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                     [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]]
                h = rad.ObjHexahedron(v, [0, 0, 0])
                rad.MatApl(h, mat_factory())
                objs.append(h)
    return objs


def get_peak_memory_mb():
    mem = psutil.Process(os.getpid()).memory_info()
    return mem.peak_wset / (1024 * 1024) if hasattr(mem, "peak_wset") else mem.rss / (1024 * 1024)


def solve_case(nside, H0, material, prec, max_iter, method):
    """material: 'bh' (MatSatIsoTab) or 'linear' (MatLin at the BH initial slope)."""
    rad.UtiDelAll()
    rad.set_demag_backend("collocation_mmmm")
    if material == "bh":
        mat_factory = lambda: rad.MatSatIsoTab(BH_DATA)
    else:
        mat_factory = lambda: rad.MatLin(MU_R_LINEAR_REF)
    objs = build_ctype(nside, mat_factory)
    cont = rad.ObjCnt(objs + [rad.ObjBckg(lambda p: [0.0, MU0 * H0, 0.0])])
    t0 = time.perf_counter()
    rad.Solve(cont, prec, max_iter, method)
    t_wall = time.perf_counter() - t0
    st = dict(rad.GetSolveStats())
    Ms = np.asarray([rad.ObjM(o)["magnetization"] for o in objs], float)
    maxM = float(np.max(np.linalg.norm(Ms, axis=1)))
    B_gap = [float(b) for b in rad.Fld(cont, "b", GAP_PROBE)]
    B_ext = [[float(b) for b in rad.Fld(cont, "b", p)] for p in EXT_PROBES]
    n_hex = len(objs)
    rad.UtiDelAll()
    return dict(
        nside=nside, n_hex=n_hex, ndof=6 * n_hex, material=material, H0_A_per_m=H0,
        method=method, prec=prec,
        t_setup=float(st.get("t_moment_system_build", 0.0)),
        t_solve=float(st.get("t_linear_solve", 0.0)),
        t_wall=t_wall,
        iterations=int(st.get("linear_iterations", -1)),
        nonl_iterations=int(st.get("nonl_iterations", -1)),
        converged=True,
        max_abs_M=maxM,
        B_gap_T=B_gap, B_gap_mag_T=float(np.linalg.norm(B_gap)),
        B_ext_T=B_ext,
        peak_memory_mb=get_peak_memory_mb(),
        num_threads=int(st.get("num_threads", 0)),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--nside", type=int, default=54, help="54 -> 27,720 hex = 166,320 DoF (~165k)")
    ap.add_argument("--calib-nside", type=int, default=24)
    ap.add_argument("--h0", type=float, default=0.0, help="drive H0 [A/m]; 0 = self-calibrate to the BH knee")
    ap.add_argument("--calib-levels", type=float, nargs="*",
                    default=[2000.0, 5000.0, 10000.0, 20000.0, 40000.0])
    ap.add_argument("--prec", type=float, default=1e-3,
                    help="nonlinear tolerance on magnetization (max|dB|/B_sat per Picard step); "
                         "1e-3 = the lab engineering standard (yano-era practice, Picard suffices)")
    ap.add_argument("--max-iter", type=int, default=300)
    ap.add_argument("--method", type=int, default=2, choices=(0, 1, 2))
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args()

    runs = {"calibration": []}
    H0 = a.h0
    if H0 <= 0.0:
        for lev in a.calib_levels:
            c = solve_case(a.calib_nside, lev, "bh", a.prec, a.max_iter, a.method)
            runs["calibration"].append(c)
            print(json.dumps(dict(calib=dict(H0=lev, max_abs_M=round(c["max_abs_M"], 1),
                                             B_gap_mag_T=round(c["B_gap_mag_T"], 5),
                                             nonl_iters=c["nonl_iterations"],
                                             wall=round(c["t_wall"], 2)))), flush=True)
            if c["max_abs_M"] >= M_KNEE:
                H0 = lev
                break
        if H0 <= 0.0:
            H0 = a.calib_levels[-1]
            print(json.dumps(dict(note=f"knee not reached; using max level {H0}")), flush=True)

    main_nl = solve_case(a.nside, H0, "bh", a.prec, a.max_iter, a.method)
    print(json.dumps(dict(main_nonlinear=main_nl)), flush=True)
    main_lin = solve_case(a.nside, H0, "linear", a.prec, a.max_iter, a.method)
    print(json.dumps(dict(linear_reference=main_lin)), flush=True)

    data = dict(
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        benchmark="moment_ctype_nonlinear_gap_field",
        problem=dict(nside=a.nside, ndof=main_nl["ndof"], BH_DATA=BH_DATA,
                     mu_r_linear_ref=MU_R_LINEAR_REF, H0_A_per_m=H0, drive_axis="y",
                     gap_probe=GAP_PROBE, method=a.method, prec=a.prec,
                     m_knee_A_per_m=M_KNEE,
                     radia_version=getattr(rad, "__version__", "unknown"),
                     python_version=platform.python_version(),
                     solver_config={k: v for k, v in rad.GetSolverConfig().items()
                                    if isinstance(v, (int, float, bool, str))}),
        results=dict(calibration=runs["calibration"], main_nonlinear=main_nl,
                     linear_reference=main_lin),
    )
    out = a.json_out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "results_moment_ctype_nonlinear.json")
    with open(out, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Results saved to {out}", flush=True)

    sat = main_nl["B_gap_mag_T"] / max(main_lin["B_gap_mag_T"], 1e-30)
    print(f"\nSUMMARY  ndof={main_nl['ndof']}  H0={H0:.0f} A/m"
          f"\n  nonlinear: |B_gap|={main_nl['B_gap_mag_T']:.4f} T  (B={main_nl['B_gap_T']})"
          f"  nonl_iters={main_nl['nonl_iterations']}  wall={main_nl['t_wall']:.1f} s"
          f"\n  linear ref: |B_gap|={main_lin['B_gap_mag_T']:.4f} T  ratio nl/lin={sat:.3f}"
          f"\n  max|M| nonlinear = {main_nl['max_abs_M']:.3e} A/m", flush=True)


if __name__ == "__main__":
    main()
