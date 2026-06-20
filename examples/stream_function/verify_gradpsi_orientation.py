"""verify_gradpsi_orientation.py -- the single-wire l>=2 fix: orient contour
loops by the consistent grad-psi winding, NOT per-loop field alignment.

bench_low_turn_geometry.py showed a Z2 (l=2 shim) coil fails as a single common
series current (~70%) at every former geometry, even though the continuous
design is ~0.9% and INDEPENDENT currents on the same contour loops reach ~0.
This script root-causes AND fixes it.

ROOT CAUSE: the production equal-current sign step (_loop_field_signs) orients
each contour loop by sign(f_k . B_target) -- i.e. flip every loop so its field
helps the target.  For a MONOTONE psi (l=1 gradient) every loop already wants
the same sense, so the flip is harmless.  For a SADDLE psi (l>=2 shim) the
natural grad-psi winding has MIXED senses (the wire reverses through the saddle);
forcing them all to "help" destroys that reversal and the single common current
cancels itself.

FIX: orient each loop so its traversal direction matches the surface current
K = n_hat x grad_s(psi) (a majority vote of d . K over the loop's segments),
then drive every loop with ONE common +current.  This is the physical
stream-function winding; sign(f.B) is recovered as the l=1 special case.

VERIFIED (this script, L=1.0 cylinder former a=0.15, DSV r=0.05, --confine abe,
nlevels=20, --eval-max 40), single COMMON series current residual (relative rms):
  target      sign(f.B)   grad-psi    independent(basis)
  Z2 (l=2)      ~0.88       ~0.056         ~0
  Gx (l=1)      ~0.010      ~0.010         ~0
grad-psi orientation fixes Z2 single-wire (88% -> 5.6%, 15x) by reversing ~33/45
loops; for Gx it is identical to sign(f.B).  So it is a strict generalisation.
This is a research probe (imports calc_streamfunction internals, no production
change); productionising it = teach _loop_field_signs / the single-stroke chain
the grad-psi orientation for multi-region psi.

Run:  python verify_gradpsi_orientation.py
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src", "radia", "panels"))
sys.path.insert(0, os.path.join(REPO, "src", "radia"))

_GEN = r"""
import os, sys
from netgen.occ import Cylinder, Sphere, Pnt, Z, OCCGeometry
from ngsolve import TaskManager
gd = sys.argv[1]
with TaskManager():
    cyl = Cylinder(Pnt(0,0,-0.5), Z, r=0.15, h=1.0)
    lat = max(cyl.faces, key=lambda f: f.mass)
    OCCGeometry(lat).GenerateMesh(maxh=0.04).Save(os.path.join(gd, "coil.vol"))
    OCCGeometry(Sphere(Pnt(0,0,0), 0.05)).GenerateMesh(maxh=0.025).Save(
        os.path.join(gd, "eval.vol"))
print("ok")
"""


def _tri_grad_K(verts, psi_v, tris):
    """Per-triangle centroid + surface-current direction K = n_hat x grad_s psi."""
    cen, K = [], []
    for (a, b, c) in tris:
        pa = verts[a]; e1 = verts[b] - pa; e2 = verts[c] - pa
        nrm = np.cross(e1, e2); A2 = np.linalg.norm(nrm)
        if A2 < 1e-30 or np.linalg.norm(e1) < 1e-30:
            continue
        t1 = e1 / np.linalg.norm(e1); nhat = nrm / A2; t2 = np.cross(nhat, t1)
        M = np.array([[e1 @ t1, e1 @ t2], [e2 @ t1, e2 @ t2]])
        rhs = np.array([psi_v[b] - psi_v[a], psi_v[c] - psi_v[a]])
        try:
            ab = np.linalg.solve(M, rhs)
        except np.linalg.LinAlgError:
            continue
        g = ab[0] * t1 + ab[1] * t2
        cen.append((pa + verts[b] + verts[c]) / 3.0)
        K.append(np.cross(nhat, g))
    return np.array(cen), np.array(K)


def run_target(cs, coil, evalv, name, target_cf, nlevels=20):
    from ngsolve import TaskManager
    ap = cs.build_argparser()
    args = ap.parse_args(["--coil-vol", coil, "--eval-vol", evalv, "--order", "1",
                          "--target-cf", target_cf, "--method", "manufacture",
                          "--eval-max", "40", "--confine", "abe",
                          "--nlevels", str(nlevels)])
    with TaskManager():
        P = cs._build_problem(args)
        psi_f, _, homo = cs._solve_and_metrics(P, args.alpha)
        psi = P["R"] @ psi_f
        coilm = P["coil"]
        verts = np.array([list(v.point) for v in coilm.vertices])
        psi_v = psi[:coilm.nv]
        tris = cs._surface_triangles(coilm)
        mpts, vb, Bm = P["mpts"], P["vector_b"], P["Bm"]

    tol = max(1e-12, 1e-6 * cs._char_len(verts, tris))
    loops = []
    for lv in cs._contour_levels(psi_v, nlevels):
        loops.extend(cs._segs_to_polylines(
            cs._contour_segments(verts, psi_v, tris, lv), tol))
    loops = [p for p in loops if len(p) >= 3]
    den = float(np.linalg.norm(Bm)) + 1e-30
    fields = [cs._wire_field_flat(cs._close_loop(p), mpts, vb, 1.0) for p in loops]

    _, signsA = cs._loop_field_signs(loops, mpts, Bm, vb)
    GA = np.sum([s * f for s, f in zip(fields, signsA)], axis=0)
    rA = float(np.linalg.norm(cs._best_fit_current(GA, Bm) * GA - Bm) / den)

    tcen, tK = _tri_grad_K(verts, psi_v, tris)
    fieldsB, n_flip = [], 0
    for p in loops:
        vote = sum(float((p[i + 1] - p[i]) @ tK[int(np.argmin(
            np.sum((tcen - 0.5 * (p[i] + p[i + 1])) ** 2, axis=1)))])
            for i in range(len(p) - 1))
        q = p if vote >= 0 else p[::-1]
        n_flip += int(vote < 0)
        fieldsB.append(cs._wire_field_flat(cs._close_loop(q), mpts, vb, 1.0))
    GB = np.sum(fieldsB, axis=0)
    rB = float(np.linalg.norm(cs._best_fit_current(GB, Bm) * GB - Bm) / den)

    Fm = np.column_stack([s * f for f, s in zip(fields, signsA)])
    Ic, *_ = np.linalg.lstsq(Fm, Bm, rcond=None)
    rC = float(np.linalg.norm(Fm @ Ic - Bm) / den)

    print(f"  {name:24s} sign(f.B)={rA:.4f}  grad-psi={rB:.4f}  "
          f"indep={rC:.4f}  (reversed {n_flip}/{len(loops)})")
    return rA, rB, rC


def main():
    import calc_streamfunction as cs
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run([sys.executable, "-c", _GEN, td],
                           capture_output=True, text=True,
                           env={**os.environ, "PYTHONIOENCODING": "utf-8"},
                           timeout=600)
        if r.returncode != 0:
            raise RuntimeError(f"mesh gen failed:\n{r.stderr[-800:]}")
        coil = os.path.join(td, "coil.vol")
        evalv = os.path.join(td, "eval.vol")
        print("single COMMON series-current residual (relative rms):")
        z2 = run_target(cs, coil, evalv, "Z2 (l=2 shim)", "z*z-(x*x+y*y)/2")
        gx = run_target(cs, coil, evalv, "Gx (l=1 control)", "x")

    ok = (z2[1] < 0.15 < z2[0]            # grad-psi fixes Z2, sign(f.B) fails it
          and abs(gx[0] - gx[1]) < 0.02)  # identical for the monotone l=1 control
    print("\nRESULT:", "CONFIRMED -- grad-psi orientation fixes single-wire l>=2"
          if ok else "NOT confirmed (investigate)")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
