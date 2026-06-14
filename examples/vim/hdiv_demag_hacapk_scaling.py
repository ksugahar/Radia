"""hdiv_demag_hacapk_scaling.py -- HACApK H-matrix compression/scaling for the HDiv-VIM charge Gram.

The HDiv-VIM demag operator N = B^T G B applies the charge-charge Coulomb Gram G (1/r over the n_charge
volume-cell + boundary-face charges) as a HACApK H-matrix (_ChargeGramHMatrix).  This benchmark answers
"does HACApK actually COMPRESS / SCALE for this problem, or is it materialize-fallback-dominated like
the compact MMM/MSC case (CLAUDE.md, 87% near-field)?" -- by building the H-matrix on a sequence of
sphere meshes (n_charge growing) and reporting the ACA stats.

RESULT (sphere, eps=1e-4, leaf=32): the compression IMPROVES with N --

  n_charge   n_lowrank   H-mat MB   dense MB   compr
       322           0       0.79       0.79    1.00   (too small: one near-field block, no compression)
       795          34       4.58       4.82    0.95
      1800         540      15.31      24.72    0.62
      3560        1780      35.26      96.69    0.37

Over N 322->3560 (11x): dense memory 122x (= N^2, confirming O(N^2)) but H-matrix memory only 45x
(~ N^1.6) -- SUB-QUADRATIC, trending to O(N log N).  ACA genuinely compresses the far blocks
(n_lowrank 0->1780).  This is BETTER than the compact MMM/MSC materialize-fallback caveat: the charge
Gram is a cleaner far-field 1/r kernel.  matvec rel err stays within the ACA eps (1e-4..1e-3; tighter
eps for large N trades some compression for accuracy).

HONEST SCOPE: demonstrated to n~3560 here (build_demag forms the dense G as the test reference, which
is O(N^2) and caps N); the TREND is clear + favorable.  Larger-N (10k+) confirmation needs a charge
extraction that bypasses the dense-G reference -- the remaining M0 scalability item.
"""
import json
import os
import time

import numpy as np
import ngsolve as ng
from netgen.csg import CSGeometry, Sphere, Pnt

import radia._radia_pybind as _rp
from radia.vim import _core as tet

HERE = os.path.dirname(os.path.abspath(__file__))


def _sphere(h):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    return ng.Mesh(g.GenerateMesh(maxh=h))


def run():
    rows = []
    with ng.TaskManager():
        for h in (0.45, 0.30, 0.22, 0.17):
            d = tet.build_demag(_sphere(h), nsub=4)
            n = d["n_charge"]
            t0 = time.perf_counter()
            H = _rp._ChargeGramHMatrix(list(d["cent"].ravel()), list(d["meas"]),
                                       list(d["self_energy"]), 1e-4, 32, 2.0)
            t_build = time.perf_counter() - t0
            st = H.stats()
            rng = np.random.default_rng(0)
            q = rng.standard_normal(n)
            t0 = time.perf_counter()
            y = np.asarray(H.matvec(q.tolist()), float)
            t_mv = time.perf_counter() - t0
            relerr = float(np.linalg.norm(y - d["G"] @ q) / np.linalg.norm(d["G"] @ q))
            rows.append(dict(maxh=h, n_charge=int(n), n_dense=int(st["n_dense"]),
                             n_lowrank=int(st["n_lowrank"]), memory_mb=float(st["memory_mb"]),
                             dense_memory_mb=float(st["dense_memory_mb"]),
                             compression=float(st["compression"]), build_s=float(t_build),
                             matvec_ms=float(t_mv * 1e3), matvec_relerr=relerr))
            print(f"  n={n:>5}  n_lowrank={st['n_lowrank']:>5}  memMB={st['memory_mb']:>7.2f}  "
                  f"denseMB={st['dense_memory_mb']:>7.2f}  compr={st['compression']:.3f}  relerr={relerr:.1e}")
    return rows


if __name__ == "__main__":
    rows = run()
    out = os.path.join(HERE, "hdiv_demag_hacapk_scaling.json")
    with open(out, "w") as f:
        json.dump({"kernel": "charge-Gram 1/r (HDiv-VIM)", "body": "unit sphere",
                   "eps": 1e-4, "leaf": 32, "rows": rows}, f, indent=2)
    print("saved", out)
    if len(rows) >= 2:
        nlo, nhi = rows[0]["n_charge"], rows[-1]["n_charge"]
        print(f"  N {nlo}->{nhi} ({nhi/nlo:.0f}x): dense mem {rows[-1]['dense_memory_mb']/rows[0]['dense_memory_mb']:.0f}x "
              f"(=N^2) vs H-matrix mem {rows[-1]['memory_mb']/rows[0]['memory_mb']:.0f}x (sub-quadratic).")
