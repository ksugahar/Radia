"""bench_yano_vs_hdiv.py -- performance comparison of two VIM (volume integral method) discretizations of
the SAME soft-iron demag problem, under TaskManager.  DEBUG-scale (h not large).

NAMING (clarified 2026-06-23): MMM, MSC, and HDiv-VIM are ALL volume integral methods (VIM) -- the unknown
is the meshed body's volume magnetization M, the kernel is the free-space Laplace Green 1/(4*pi*r), and the
open boundary is natural (no air mesh / Kelvin / PML).  "VIM" is the method CLASS, not a single method.  The
members differ only in the M discretization + projection:
  * MSC (yano-type) : surface charge sigma per face, COLLOCATION (EIEM2).  hex 6 DOF / wedge 5 DOF.
  * HDiv (FEEC)     : M in H(div) RT_k, GALERKIN  (N = B^T G B, the C++ charge-Gram).  1 flux/face at RT0.
  (MMM is the moment/dipole collocation member on tets.)  So this bench is "MSC vs HDiv", two VIM members --
  NOT "MSC vs VIM" (that mixes a member with the class).

WHAT THIS EXERCISES
  * MSC (yano)  : rad.Solve(.., demag_backend="yano", solver=2) -> the C++ HACApK MSC method-2 BiCGSTAB
                  (HEX-only, 6 sigma/hex).  Self-wrapped in a RegionTaskManager (CLAUDE.md "C++ HACApK
                  Self-Wrap Policy", 2026-06-23) -> parallel build + solve.
  * HDiv (FEEC) : rad.Solve(.., demag_backend="hdiv") -> radia.vim.hdiv_demag_solve (RT0) on the SAME hex
                  iron, C++ _ChargeGramHMatrix demag.  The Python dense O(N^2) charge-Gram was REMOVED
                  2026-06-23 -- the C++ H-matrix kernel IS the demag operator, so THIS bench times C++, not
                  a "clever Python matrix build".

TWO AXES
  Part 1 (HEX, h-sweep)  : MSC vs HDiv-RT0 self-consistent demag solve.  t_total, iters, ndof, and an
                           EXTERNAL probe B_z (the agreement metric -- compare external field, never the
                           inside-material field) so the two methods are checked to compute the same physics.
  Part 2 (TET, p-sweep)  : HDiv p-convergence via radia.vim.DemagOperator(HDiv(mesh, order=p)).  The
                           high-order monomial-charge Gram is TET-only and MSC is HEX-only (6 DOF/hex), so
                           each method's convergence axis uses the geometry it supports.

Run:  python examples/vim/bench_yano_vs_hdiv.py        (full sweep; saves bench_yano_vs_hdiv.json here)
Fast: BENCH_N="2,3" BENCH_P="0,1" python examples/vim/bench_yano_vs_hdiv.py
NOTE: this loads the shared radia .pyd; if a parallel build is rebuilding it, runs can crash/hang -- run on
a quiet machine (or against a pinned .pyd snapshot) for clean timing.
"""
import json
import os
import platform
import time

import ngsolve as ng
from ngsolve.meshes import MakeStructured3DMesh
from netgen.occ import Box, OCCGeometry, Pnt

import radia as rad
import radia.vim as vim
from radia.vim import DemagOperator

HERE = os.path.dirname(os.path.abspath(__file__))
N_THREADS = int(os.environ.get("BENCH_THREADS", "4"))
ng.SetNumThreads(N_THREADS)

N_LIST = [int(x) for x in os.environ.get("BENCH_N", "2,3,4,5,6").split(",")]
P_LIST = [int(x) for x in os.environ.get("BENCH_P", "0,1,2,3").split(",")]
TET_MAXH = float(os.environ.get("BENCH_MAXH", "0.5"))

MU_R = 100.0
H0 = 1.0e3                                  # A/m (linear -> magnitude is timing-irrelevant)
PROBE = [0.5, 0.5, 2.0]                     # external probe, 1 cube-length above the top face centre


def _peak_mem_mb():
    try:
        import psutil
        m = psutil.Process(os.getpid()).memory_info()
        return (getattr(m, "peak_wset", None) or m.rss) / (1024 * 1024)
    except Exception:
        return float("nan")


def _bkg():
    return rad.ObjBckg(lambda p: [0.0, 0.0, H0])


def part1_hex_hsweep(n_list):
    print("=" * 80, flush=True)
    print("Part 1  HEX cube  --  MSC (yano-type collocation) vs HDiv (RT0, FEEC Galerkin)", flush=True)
    print(f"        mu_r={MU_R}, H0={H0} A/m, {N_THREADS} threads, external probe B_z at {PROBE}", flush=True)
    print("=" * 80, flush=True)
    rows = []
    for n in n_list:
        mesh = MakeStructured3DMesh(hexes=True, nx=n, ny=n, nz=n)
        nhex = mesh.GetNE(ng.VOL)

        rad.UtiDelAll()
        iron = vim.soft_iron_from_mesh(mesh, mu_r=MU_R)
        model = rad.ObjCnt([iron, _bkg()])
        with ng.TaskManager():
            t0 = time.perf_counter()
            rad.Solve(model, 1e-5, 2000, 2, demag_backend="yano")
            t_msc = time.perf_counter() - t0
            Bz_msc = float(rad.Fld(model, "b", PROBE)[2])
        ndof_msc = 6 * nhex

        rad.UtiDelAll()
        iron = vim.soft_iron_from_mesh(mesh, mu_r=MU_R)
        model = rad.ObjCnt([iron, _bkg()])
        with ng.TaskManager():
            t0 = time.perf_counter()
            res = rad.Solve(model, 1e-5, 2000, 2, demag_backend="hdiv")
            t_hdiv = time.perf_counter() - t0
            Bz_hdiv = float(rad.Fld(model, "b", PROBE)[2])

        dBz = 100.0 * (Bz_hdiv - Bz_msc) / Bz_msc if Bz_msc else float("nan")
        print(f"  n={n}^3 nhex={nhex:<4} | MSC : ndof={ndof_msc:<5} t={t_msc:7.3f}s Bz={Bz_msc:.5e} "
              f"| HDiv: ndof={res['ndof']:<5} t={t_hdiv:7.3f}s iters={res['iters']:<3} "
              f"demag={res['demag']:.4f} Bz={Bz_hdiv:.5e} | dBz={dBz:+.3f}%", flush=True)
        rows.append(dict(n=n, nhex=int(nhex),
                         msc=dict(ndof=ndof_msc, t_total=t_msc, Bz_ext=Bz_msc),
                         hdiv=dict(ndof=int(res["ndof"]), t_total=t_hdiv, iters=int(res["iters"]),
                                   demag=float(res["demag"]), M_avg_z=float(res["M_avg"][2]), Bz_ext=Bz_hdiv),
                         dBz_ext_pct=dBz, speedup_msc_over_hdiv=(t_hdiv / t_msc if t_msc else float("nan"))))
    return rows


def part2_tet_psweep(maxh, p_list):
    print("\n" + "=" * 80, flush=True)
    print(f"Part 2  TET cube (maxh={maxh})  --  HDiv (FEEC) p-convergence: demag factor -> ~1/3", flush=True)
    print("=" * 80, flush=True)
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=maxh))
    ntet = mesh.GetNE(ng.VOL)
    print(f"  ntet={ntet}", flush=True)
    rows = []
    Mz = ng.CoefficientFunction((0.0, 0.0, 1.0))
    for p in p_list:
        with ng.TaskManager():
            fes = ng.HDiv(mesh, order=p)
            t0 = time.perf_counter()
            op = DemagOperator(fes, eps=1e-7)               # BUILD the C++ charge-Gram H-matrix
            t_build = time.perf_counter() - t0
            t0 = time.perf_counter()
            D = op.DemagFactor(Mz)
            t_factor = time.perf_counter() - t0
            ndof = op.ndof
        print(f"  p={p} ndof={ndof:<5} t_build={t_build:7.3f}s t_factor={t_factor:.4f}s demag={D:.5f}", flush=True)
        rows.append(dict(p=p, ndof=int(ndof), t_build=t_build, t_factor=t_factor, demag=float(D)))
    return dict(maxh=maxh, ntet=int(ntet), orders=rows)


def main():
    t_start = time.perf_counter()
    part1 = part1_hex_hsweep(N_LIST)
    part2 = part2_tet_psweep(TET_MAXH, P_LIST)
    out = dict(benchmark="msc_yano_vs_hdiv_feec", hostname=platform.node(), n_threads=N_THREADS,
               problem=dict(geometry="unit iron cube", mu_r=MU_R, H0_Am=H0, probe=PROBE),
               peak_memory_mb=_peak_mem_mb(), wall_time_s=time.perf_counter() - t_start,
               part1_hex_hsweep=part1, part2_tet_psweep=part2,
               notes=("Both are VIM members: MSC=yano-type collocation (HEX), HDiv=FEEC Galerkin.  C++ "
                      "_ChargeGramHMatrix is the sole demag operator (Python dense Gram removed 2026-06-23). "
                      "All solves under TaskManager (HACApK build+solve self-wrap)."))
    path = os.path.join(HERE, "bench_yano_vs_hdiv.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"\npeak_mem={out['peak_memory_mb']:.0f} MB, wall={out['wall_time_s']:.1f}s\nsaved {path}", flush=True)


if __name__ == "__main__":
    main()
