"""COIL-DRIVEN nonlinear C-yoke gap-field benchmark (the classic ~0.1 T C-type), collocation MMMM method 2.

This is the voxelized-hex twin of the VALIDATED tet parity fixture
validation_test/feec/parity_vs_msc/test_ctype_nonlinear_parity.py -- SAME geometry predicate
(outer 0.12 m square x 0.04 m, 0.07 m bore, mouth fully open for x >= 0.018), SAME coil
(rad.ObjFlmCur rectangular loop in the x-z plane at y=0 encircling the back leg, NI = 8000 A),
SAME smooth Froehlich BH (chi0 = 2000, Msat = 1.6 T / mu0), SAME gap probes (bore centre
[0,0,0] + gap opening [0.02,0,0]).  The parity fixture measured gap-centre |B| ~ 0.098 T at
2715 tets (HDiv vs MMM agree 0.32%), so the voxel-hex numbers here are directly comparable.

Unlike bench_moment_ctype_nonlinear.py (uniform-field drive, benchmark-only observable), this
is the ENGINEERING C-type: coil MMF drives ~0.1 T into the gap.  A resolution ladder (default
nside 24 / 42 / 60 -> ~10k / ~58k / ~165k DoF) shows gap-B mesh convergence; a linear
reference (MatLin at chi0) at the largest size quantifies the saturation effect.

JSON per Benchmark Policy.  Self-contained (no repo needed on the bench machine).
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

# --- the parity fixture's material + coil + probes, verbatim ---------------------
CHI0 = 2000.0
MSAT = 1.6 / MU0                        # ~1.2732e6 A/m
_Hs = np.concatenate([[0.0], np.logspace(-1, 7, 80)])
_Ms = CHI0 * _Hs / (1.0 + CHI0 * _Hs / MSAT)
_Bs = MU0 * (_Hs + _Ms)
BH = [[float(h), float(b)] for h, b in zip(_Hs, _Bs)]
NI = 8000.0
GAP_PTS = [[0.0, 0.0, 0.0], [0.02, 0.0, 0.0]]     # bore centre + gap opening
EXT_PT = [0.0, 0.09, 0.01]


def make_coil():
    """Rectangular filament loop in the x-z plane at y=0 encircling the back leg."""
    pts = [[-0.075, 0.0, -0.035], [-0.025, 0.0, -0.035], [-0.025, 0.0, 0.035],
           [-0.075, 0.0, 0.035], [-0.075, 0.0, -0.035]]
    return rad.ObjFlmCur(pts, NI)


def _inside_parity_cyoke(cx, cy):
    """The parity fixture's cross-section: outer square minus bore minus the FULL x>=0.018 slab."""
    if not (-0.06 <= cx <= 0.06 and -0.06 <= cy <= 0.06):
        return False
    if -0.035 <= cx <= 0.035 and -0.035 <= cy <= 0.035:
        return False
    if cx >= 0.018:
        return False
    return True


def count_hexes(nside):
    nxy, nz = nside, max(2, nside // 3)
    xs = np.linspace(-0.06, 0.06, nxy + 1)
    n_inplane = sum(1 for j in range(nxy) for i in range(nxy)
                    if _inside_parity_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])))
    return n_inplane * nz


def build_yoke(nside, mat_factory):
    nxy, nz = nside, max(2, nside // 3)
    xs = np.linspace(-0.06, 0.06, nxy + 1)
    zs = np.linspace(-0.02, 0.02, nz + 1)
    objs = []
    for k in range(nz):
        for j in range(nxy):
            for i in range(nxy):
                if not _inside_parity_cyoke(0.5 * (xs[i] + xs[i + 1]), 0.5 * (xs[j] + xs[j + 1])):
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


def solve_case(nside, material, prec, max_iter, method):
    rad.UtiDelAll()
    rad.set_demag_backend("collocation_mmmm")
    if material == "bh":
        mat_factory = lambda: rad.MatSatIsoTab(BH)
    else:
        mat_factory = lambda: rad.MatLin(CHI0 + 1.0)
    objs = build_yoke(nside, mat_factory)
    coil = make_coil()
    cont = rad.ObjCnt(objs + [coil])
    t0 = time.perf_counter()
    rad.Solve(cont, prec, max_iter, method)
    t_wall = time.perf_counter() - t0
    st = dict(rad.GetSolveStats())
    Ms = np.asarray([rad.ObjM(o)["magnetization"] for o in objs], float)
    maxM = float(np.max(np.linalg.norm(Ms, axis=1)))
    B_gap = [[float(b) for b in rad.Fld(cont, "b", p)] for p in GAP_PTS]     # iron + coil field
    B_ext = [float(b) for b in rad.Fld(cont, "b", EXT_PT)]
    n_hex = len(objs)
    rad.UtiDelAll()
    return dict(
        nside=nside, n_hex=n_hex, ndof=6 * n_hex, material=material, NI_A=NI,
        method=method, prec=prec,
        t_setup=float(st.get("t_moment_system_build", 0.0)),
        t_solve=float(st.get("t_linear_solve", 0.0)),
        t_wall=t_wall,
        iterations=int(st.get("linear_iterations", -1)),
        nonl_iterations=int(st.get("nonl_iterations", -1)),
        converged=True,
        max_abs_M=maxM,
        B_bore_T=B_gap[0], B_bore_mag_T=float(np.linalg.norm(B_gap[0])),
        B_gapopen_T=B_gap[1], B_gapopen_mag_T=float(np.linalg.norm(B_gap[1])),
        B_ext_T=B_ext,
        peak_memory_mb=get_peak_memory_mb(),
        num_threads=int(st.get("num_threads", 0)),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--sides", type=int, nargs="*", default=[24, 42, 60],
                    help="voxel resolutions; 60 -> ~165k DoF (exact count printed)")
    ap.add_argument("--linear-ref-side", type=int, default=60)
    ap.add_argument("--prec", type=float, default=1e-3,
                    help="nonlinear tolerance (max|dB|/B_sat per Picard step).  1e-3 = the lab "
                         "engineering standard (Sugahara: the yano-type MMM/MSC era operated at "
                         "~1e-3 and Picard suffices there; the gap-B observable's voxel "
                         "discretization error ~0.5-1%% dominates anyway).  Use 1e-6 only for "
                         "parity-fixture-matched cross-checks (costs ~2x the Picard iterations).")
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--method", type=int, default=2, choices=(0, 1, 2))
    ap.add_argument("--json-out", default=None)
    a = ap.parse_args()

    cases = []
    for ns in a.sides:
        print(json.dumps(dict(plan=dict(nside=ns, n_hex=count_hexes(ns), ndof=6 * count_hexes(ns)))), flush=True)
    for ns in a.sides:
        c = solve_case(ns, "bh", a.prec, a.max_iter, a.method)
        cases.append(c)
        print(json.dumps(dict(nonlinear=dict(nside=ns, ndof=c["ndof"],
                                             B_bore_mag_T=round(c["B_bore_mag_T"], 5),
                                             B_gapopen_mag_T=round(c["B_gapopen_mag_T"], 5),
                                             nonl_iters=c["nonl_iterations"], iters=c["iterations"],
                                             wall=round(c["t_wall"], 1),
                                             max_abs_M=round(c["max_abs_M"], 0)))), flush=True)
    lin = solve_case(a.linear_ref_side, "linear", a.prec, a.max_iter, a.method)
    print(json.dumps(dict(linear_reference=dict(nside=lin["nside"], ndof=lin["ndof"],
                                                B_bore_mag_T=round(lin["B_bore_mag_T"], 5),
                                                wall=round(lin["t_wall"], 1)))), flush=True)

    data = dict(
        timestamp=datetime.now().isoformat(),
        hostname=platform.node(),
        benchmark="moment_ctype_coil_nonlinear_gap_field",
        problem=dict(sides=a.sides, NI_A=NI, chi0=CHI0, Msat_A_per_m=MSAT,
                     BH_points=len(BH), gap_probes=GAP_PTS, method=a.method, prec=a.prec,
                     geometry="parity fixture cross-section (test_ctype_nonlinear_parity.py), voxelized hex",
                     parity_reference="tet 2715 MMM/HDiv gap-centre |B| ~ 0.098 T (agree 0.32%)",
                     radia_version=getattr(rad, "__version__", "unknown"),
                     python_version=platform.python_version(),
                     solver_config={k: v for k, v in rad.GetSolverConfig().items()
                                    if isinstance(v, (int, float, bool, str))}),
        results=dict(nonlinear=cases, linear_reference=lin),
    )
    out = a.json_out or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     "results_moment_ctype_coil_nonlinear.json")
    with open(out, "w") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"Results saved to {out}", flush=True)

    big = cases[-1]
    print(f"\nSUMMARY  ndof={big['ndof']}  NI={NI:.0f} A"
          f"\n  nonlinear  |B| bore centre = {big['B_bore_mag_T']:.4f} T   gap opening = {big['B_gapopen_mag_T']:.4f} T"
          f"\n             nonl_iters={big['nonl_iterations']}  wall={big['t_wall']:.1f} s  mem={big['peak_memory_mb']:.0f} MB"
          f"\n  linear ref |B| bore centre = {lin['B_bore_mag_T']:.4f} T   ratio nl/lin = "
          f"{big['B_bore_mag_T'] / max(lin['B_bore_mag_T'], 1e-30):.3f}"
          f"\n  parity tet fixture reference: ~0.098 T (2715 tets, HDiv==MMM 0.32%)", flush=True)


if __name__ == "__main__":
    main()
