"""Cohomology generators FROM SCRATCH via the matched Whitney complex.

The net-current (secular) basis of a stream-function coil designer is the first
cohomology H^1 of the winding surface (dim = first Betti number b1).  Instead of
hard-coding analytic angle-gradients n x grad(theta) -- which are valid div-free
class representatives but NOT the canonical harmonic forms (curl_s(n x grad zeta)
= Delta_s zeta != 0 on a curved surface) -- this computes the generators
NUMERICALLY, with no analytic ansatz, from the EXACT discrete de Rham complex:

    H1(order=1, P1 nodal)  --grad-->  HCurl(order=0, Whitney edges)  --curl-->  L2(P0)

The Whitney pairing makes curl(grad)=0 to machine precision (the MATCHED complex;
a mismatched HCurl(1)/H1(1) pairing breaks exactness and the harmonics cannot
separate).  The harmonic 1-forms are then the kernel of the Hodge-Laplacian proxy

    K = curl-curl  +  grad-div-penalty  =  Acc + B M0^{-1} B^T

(both terms vanish only on curl-free AND div-free fields), found as the b1
smallest generalised eigenpairs K x = lam M x.  They are the canonical
(minimum-norm) secular currents -- the lowest-complexity net-current coils.

Validation: (1) the complex is exact (curl of a gradient ~ 0); (2) exactly b1
near-zero eigenvalues with a huge gap; (3) each generator is div-free AND
curl-free (harmonic); (4) they reproduce the analytic generators' NET-CURRENT
(TF Ampere) field -- the SAME cohomology class; (5) the surface Euler
characteristic (gmsh-free) confirms b1.

Run:  python validation_test/stream_function/cohomology_generators_whitney.py
"""
from __future__ import annotations

import math
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
_EXAMPLES = os.path.join(_REPO, "examples", "stream_function")
sys.path.insert(0, _HERE)
sys.path.insert(0, _EXAMPLES)
sys.path.insert(0, os.path.join(_REPO, "src", "radia"))
sys.path.insert(0, os.path.join(_REPO, "src", "radia", "panels"))


def _to_scipy(mat):
    r, c, v = mat.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))))


def main(no_cohomology=False):
    import demo_regcoil_fusion as D
    from ngsolve import (Mesh, HCurl, H1, specialcf, TaskManager, BilinearForm,
                         curl, grad, ds, InnerProduct, GridFunction, Integrate,
                         x, y, z, sqrt)

    print("=== cohomology generators from scratch (matched Whitney complex) ===")
    work = os.path.join(os.environ.get("TEMP", "/tmp"), "sf_coh_whitney")
    os.makedirs(work, exist_ok=True)

    with TaskManager():
        coil = Mesh(D._torus_surface_vol(D.A_WIND, 0.05,
                                         os.path.join(work, "wind.vol")))
        n_cf = specialcf.normal(3)
        # MATCHED Whitney pair: HCurl order 0 (edges), H1 order 1 (nodes)
        fes_h = HCurl(coil, order=0, definedon=coil.Boundaries(".*"))
        fes_0 = H1(coil, order=1, definedon=coil.Boundaries(".*"))
        u, v = fes_h.TnT()
        phi, q = fes_0.TnT()

        acc = BilinearForm(fes_h)
        acc += curl(u).Trace() * curl(v).Trace() * ds
        acc.Assemble()
        mm = BilinearForm(fes_h)
        mm += InnerProduct(u.Trace(), v.Trace()) * ds
        mm.Assemble()
        bf = BilinearForm(trialspace=fes_0, testspace=fes_h)
        bf += InnerProduct(v.Trace(), grad(phi).Trace()) * ds
        bf.Assemble()
        m0 = BilinearForm(fes_0)
        m0 += phi * q * ds
        m0.Assemble()

        Acc = _to_scipy(acc.mat); M = _to_scipy(mm.mat).tocsc()
        B = _to_scipy(bf.mat); M0 = _to_scipy(m0.mat)
        print(f"Whitney HCurl(0) ndof={fes_h.ndof} (edges), "
              f"H1(1) ndof={fes_0.ndof} (nodes)")

        # (1) exactness of the complex: curl-energy of a gradient ~ 0
        rgrad = (np.arange(fes_0.ndof) % 7 - 3.0)
        gfield = spla.spsolve(M, B @ rgrad)
        exact_resid = abs(float(gfield @ (Acc @ gfield))) / (
            float(gfield @ (M @ gfield)) + 1e-30)

        # Hodge-Laplacian proxy kernel = harmonic forms
        d0 = np.maximum(np.array(M0.sum(axis=1)).ravel(), 1e-30)
        P = B @ sp.diags(1.0 / d0) @ B.T
        K = (Acc + P).tocsc()
        vals, vecs = spla.eigsh(K + 1e-12 * M, k=6, M=M, sigma=0.0, which="LM")
        o = np.argsort(np.abs(vals)); vals = vals[o]; vecs = vecs[:, o]
        gap = abs(vals[2]) / (abs(vals[1]) + 1e-30)
        b1_numeric = int(np.sum(np.abs(vals) < 1e-6 * abs(vals[2])))

        # (3) each generator: div-free + curl-free (Rayleigh quotients ~ 0)
        gens = []
        gen_checks = []
        for j in range(2):
            xj = vecs[:, j].real
            mj = float(xj @ (M @ xj))
            curl_e = float(xj @ (Acc @ xj)) / (mj + 1e-30)
            div_e = float(xj @ (P @ xj)) / (mj + 1e-30)
            gen_checks.append({"curl_e": curl_e, "div_e": div_e})
            gf = GridFunction(fes_h)
            gf.vec.FV().NumPy()[:] = xj
            gens.append(gf)

        # (4) class match: fit analytic generators' NET-CURRENT (TF Ampere) field
        def tf_BR(Kc):
            radii = np.linspace(0.24, 0.36, 4)
            by = np.zeros(len(radii))
            for i, r in enumerate(radii):
                dxt, dyt, dzt = r - x, -y, -z
                r3 = (dxt*dxt+dyt*dyt+dzt*dzt) * sqrt(dxt*dxt+dyt*dyt+dzt*dzt)
                by[i] = 1e-7 * Integrate((Kc[2]*dxt - Kc[0]*dzt) / r3 * ds, coil)
            return by * radii
        tf_analytic = tf_BR(D._secular_currents(n_cf)[0][1])
        Hfields = np.column_stack([tf_BR(g) for g in gens])
        coef, *_ = np.linalg.lstsq(Hfields, tf_analytic, rcond=None)
        class_resid = float(np.linalg.norm(Hfields @ coef - tf_analytic)
                            / (np.linalg.norm(tf_analytic) + 1e-30))

        b1_euler = None if no_cohomology else D._betti1_winding_surface(D.A_WIND,
                                                                        0.06)

    print(f"[exact ] curl-energy of a gradient = {exact_resid:.2e}  (~0)")
    print(f"[eigs  ] lam = {', '.join(f'{v:.2e}' for v in vals[:4])} ...")
    print(f"[eigs  ] gap |lam2|/|lam1| = {gap:.2e}  -> numeric b1 = {b1_numeric}")
    for j, c in enumerate(gen_checks):
        print(f"[gen {j} ] curl-energy={c['curl_e']:.2e}  div-energy={c['div_e']:.2e}"
              f"  (both ~0 = harmonic)")
    print(f"[class ] analytic TF net-current fit by numeric gens = {class_resid:.2e}")
    print(f"[euler ] b1 = {b1_euler}  (expect 2)")

    result = {
        "n_edges": int(fes_h.ndof), "n_nodes": int(fes_0.ndof),
        "exactness_resid": exact_resid,
        "eig_gap": float(gap), "b1_numeric": b1_numeric, "b1_euler": b1_euler,
        "gen_curl_energy": [c["curl_e"] for c in gen_checks],
        "gen_div_energy": [c["div_e"] for c in gen_checks],
        "class_match_resid": class_resid,
    }
    ok = (exact_resid < 1e-8
          and b1_numeric == 2 and gap > 1e6
          and all(c["curl_e"] < 1e-6 and c["div_e"] < 1e-6 for c in gen_checks)
          and class_resid < 0.05
          and (no_cohomology or b1_euler == 2))
    print(f"\n=== PASS: {ok} ===")
    return result, ok


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-cohomology", action="store_true")
    a = ap.parse_args()
    _, ok = main(no_cohomology=a.no_cohomology)
    raise SystemExit(0 if ok else 1)
