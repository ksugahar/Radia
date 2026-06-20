# Benchmark: free-space Biot-Savart, three ways, equal target accuracy ~1e-4, scaling N=Q = 500..8000.
#
# This is the empirical backing for the CLAUDE.md "WAIT for Hlib in ngsolve.bem" policy: it measures
# where FMM (the backend ngsolve.bem ships TODAY) wins vs where an H-matrix (the Hlib/ACA backend
# expected to land in ngsolve.bem) wins, on the SF-coil free-space Biot-Savart source.
#
# Geometry = the SF-coil regime: N tangential filament segments on r=0.40 (the psi=S2 coil) -> Q field
# points on r=0.70 (the iron shell).  Concentric, radially separated (FMM-admissible by multipole
# radius; geometrically NON-separated bounding boxes -> the worst case for a single-block ACA, which
# is why the rank is ~900 -- a hierarchical HACApK would split directionally and lower per-block rank).
#
#   direct : dense numpy Biot-Savart, O(N*Q) time, O(N*Q) memory
#   FMM    : ngsolve.bem.BiotSavartCF (multipole), build O(N) + vectorised batch eval at Q points
#   ACA+   : radia.stream_function.aca_tsvd on the (3Q)x(3N) Biot-Savart block, matvec O((N+Q)*rank)
#
# Findings (LAB, see act8_14_bench_fmm_vs_aca_biotsavart.json / .png next to this file):
#   * direct O(N^2); FMM total O(N) (one-shot winner); ACA+ rank ~900 constant in N -> O(N*r) apply+mem.
#   * one-shot field eval -> FMM (low setup);  repeated apply -> ACA+/H-matrix (near-free matvec).
#   * CAVEAT: the ACA+ BUILD time here is inflated by the Python entry-callback (~1.9e7 calls); a C++
#     HACApK assembly would be far faster.  The apply / rank / memory numbers ARE algorithm-representative.
import os, sys, time, json, platform
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
from datetime import datetime
from math import sqrt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../src/radia"))
import ngsolve as ng
from ngsolve import TaskManager
import netgen.occ as occ
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve.bem import BiotSavartCF
from ngsolve.bla import Vec3D
from stream_function import aca_tsvd

INV4PI = 1.0 / (4.0 * np.pi)
MP_ORDER, KAPPA, MP_RAD = 24, 1e-6, 0.42       # FMM: order for ~1e-4 at r=0.70 ((0.42/0.70)^24 ~ 5e-6)
ACA_EPS = 1e-4
R_SRC, R_TGT = 0.40, 0.70
SIZES = [500, 1000, 2000, 4000, 8000]


def fib_sphere(n, R):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); gold = np.pi * (1 + 5 ** 0.5); th = gold * i
    return R * np.stack([np.sin(phi) * np.cos(th), np.sin(phi) * np.sin(th), np.cos(phi)], 1)


def setup(N, Q):
    C = fib_sphere(N, R_SRC)
    t = np.cross(C, np.array([0, 0, 1.0]))
    t /= np.linalg.norm(t, axis=1, keepdims=True) + 1e-30      # tangential (coil-like) direction
    L = 0.01
    SP, EP = C - 0.5 * L * t, C + 0.5 * L * t
    MOM = (EP - SP)                                            # I=1 current-element vectors
    Ce = 0.5 * (SP + EP)
    T = fib_sphere(Q, R_TGT)
    return SP, EP, Ce, MOM, T


def direct(Ce, MOM, T):
    out = np.zeros((len(T), 3))
    for i, p in enumerate(T):
        d = p - Ce; r3 = np.einsum('ij,ij->i', d, d) ** 1.5 + 1e-300
        out[i] = INV4PI * (np.cross(MOM, d) / r3[:, None]).sum(0)
    return out


def fmm(BG, SP, EP, T):
    bs = BiotSavartCF(MP_ORDER, KAPPA, Vec3D(0, 0, 0), MP_RAD)
    t0 = time.perf_counter()
    for k in range(len(SP)):
        bs.AddCurrent(Vec3D(*SP[k]), Vec3D(*EP[k]), 1.0, num=1)
    t_build = time.perf_counter() - t0
    Hcf = ng.CF((bs[0].real, bs[1].real, bs[2].real))
    t0 = time.perf_counter()
    mp = BG(T[:, 0], T[:, 1], T[:, 2])
    vals = np.array(Hcf(mp)).reshape(len(T), 3)
    t_eval = time.perf_counter() - t0
    return vals, t_build, t_eval


def aca(Ce, MOM, T):
    Q, N = len(T), len(Ce)
    SRCl = [(float(a), float(b), float(c)) for a, b, c in Ce]   # plain floats for a fast entry callback
    TGTl = [(float(a), float(b), float(c)) for a, b, c in T]

    def entry(I, J):
        i, a = I // 3, I % 3; j, b = J // 3, J % 3
        tx, ty, tz = TGTl[i]; sx, sy, sz = SRCl[j]
        dx, dy, dz = tx - sx, ty - sy, tz - sz
        r2 = dx * dx + dy * dy + dz * dz; r3 = r2 * sqrt(r2)
        # K[a,b] = -INV4PI * skew(d)[a,b] / r^3 ; skew(d) = [[0,-dz,dy],[dz,0,-dx],[-dy,dx,0]]
        if a == 0:   s = (0.0, -dz, dy)[b]
        elif a == 1: s = (dz, 0.0, -dx)[b]
        else:        s = (-dy, dx, 0.0)[b]
        return -INV4PI * s / r3

    t0 = time.perf_counter()
    res = aca_tsvd(3 * Q, 3 * N, entry, aca_eps=ACA_EPS, kmax=min(3 * Q, 3 * N))
    t_build = time.perf_counter() - t0
    mom = MOM.reshape(-1)                                      # (3N,)
    t0 = time.perf_counter()
    U, S, V = res.U[:, :res.modes], res.S[:res.modes], res.V[:, :res.modes]
    field = U @ (S * (V.T @ mom))                             # (3Q,)
    t_eval = time.perf_counter() - t0
    return field.reshape(Q, 3), t_build, t_eval, int(res.k_aca), int(res.modes)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    with TaskManager():
        box = Box(Pnt(-1.2, -1.2, -1.2), Pnt(1.2, 1.2, 1.2)); box.mat("air")
        BG = ng.Mesh(OCCGeometry(box).GenerateMesh(maxh=0.5))   # point-location host for FMM eval (untimed)

    rows = []
    for n in SIZES:
        N = Q = n
        SP, EP, Ce, MOM, T = setup(N, Q)
        t0 = time.perf_counter(); Bd = direct(Ce, MOM, T); t_dir = time.perf_counter() - t0
        Bf, tfb, tfe = fmm(BG, SP, EP, T)
        Ba, tab, tae, rank, modes = aca(Ce, MOM, T)
        ef = float(np.linalg.norm(Bf - Bd) / np.linalg.norm(Bd))
        ea = float(np.linalg.norm(Ba - Bd) / np.linalg.norm(Bd))
        dense_mb = (3 * Q) * (3 * N) * 8 / 1e6
        aca_mb = ((3 * Q) + (3 * N)) * rank * 8 / 1e6
        rows.append(dict(n=n, t_direct=t_dir,
                         fmm_build=tfb, fmm_eval=tfe, fmm_total=tfb + tfe, fmm_err=ef,
                         aca_build=tab, aca_eval=tae, aca_total=tab + tae, aca_err=ea,
                         aca_rank=rank, aca_modes=modes, dense_mb=dense_mb, aca_mb=aca_mb))
        print("N=Q=%5d | direct %7.3fs | FMM %6.3f+%6.3f=%6.3fs err %.1e | ACA %7.3f+mv %6.4f=%7.3fs "
              "rank %4d err %.1e | mem dense %8.1f aca %6.1f MB"
              % (n, t_dir, tfb, tfe, tfb + tfe, ef, tab, tae, tab + tae, rank, ea, dense_mb, aca_mb),
              flush=True)

    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="fmm_vs_aca_biotsavart",
               problem=dict(kernel="free-space Biot-Savart", R_src=R_SRC, R_tgt=R_TGT,
                            mp_order=MP_ORDER, aca_eps=ACA_EPS, sizes=SIZES),
               note="ACA+ build inflated by Python entry callback; apply/rank/memory are algorithm-representative.",
               results=rows)
    jpath = os.path.join(here, "act8_14_bench_fmm_vs_aca_biotsavart.json")
    with open(jpath, "w") as f:
        json.dump(out, f, indent=2)
    print("saved", jpath)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return
    n = np.array([r["n"] for r in rows], float)
    g = lambda k: np.array([r[k] for r in rows], float)
    fig, ax = plt.subplots(1, 2, figsize=(9.2, 3.8))
    ax[0].loglog(n, g("t_direct"), "o-", color="#444", label="direct dense  O(N^2)")
    ax[0].loglog(n, g("fmm_total"), "s-", color="#1f77b4", label="FMM total (build+eval)")
    ax[0].loglog(n, g("fmm_build"), "s--", color="#1f77b4", alpha=.5, label="FMM build")
    ax[0].loglog(n, g("aca_build"), "^-", color="#d62728", label="ACA+ build (Py-callback)")
    ax[0].loglog(n, g("aca_eval"), "^--", color="#d62728", alpha=.7, label="ACA+ apply (matvec)")
    ax[0].loglog(n, g("t_direct")[0] * (n / n[0]) ** 2, ":", color="#aaa", lw=1)
    ax[0].loglog(n, g("fmm_total")[0] * (n / n[0]) ** 1, ":", color="#9cf", lw=1)
    ax[0].set_xlabel("N = Q  (sources = field points)"); ax[0].set_ylabel("wall time  [s]")
    ax[0].legend(fontsize=7, loc="upper left"); ax[0].grid(True, which="both", alpha=.25)
    ax[1].loglog(n, g("dense_mb"), "o-", color="#444", label="dense matrix  O(N^2)")
    ax[1].loglog(n, g("aca_mb"), "^-", color="#d62728", label="ACA+ factors  O(N*r), r~900")
    ax[1].set_xlabel("N = Q"); ax[1].set_ylabel("storage  [MB]")
    ax[1].legend(fontsize=8, loc="upper left"); ax[1].grid(True, which="both", alpha=.25)
    fig.tight_layout()
    ppath = os.path.join(here, "act8_14_bench_fmm_vs_aca_biotsavart.png")
    fig.savefig(ppath, dpi=130, bbox_inches="tight")
    print("saved", ppath)


if __name__ == "__main__":
    main()
