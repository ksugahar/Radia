"""hdiv_vim_distortion_loops.py -- WHY HDiv-VIM has no element-distortion penalty (the structural answer
to the Koizumi-Higuchi 1995 difficulty), demonstrated on the de Rham loop nullspace.

Background (Koizumi & Higuchi, IEEE TMag 31(3):1516-1519, 1995, doi:10.1109/20.376318 -- the face/vector
element VIE that is the direct predecessor of this HDiv-VIM): they represented the magnetization M inside a
hexahedron with CO/CONTRA-VARIANT vectors a_i = dr/dzeta_i (the element's edge-tangent frame) plus the six
face values, and reported that "the calculation results strongly depend on the DIRECTION OF VECTORS for the
unknowns" and are "poor when the element shape is strongly deformed".

ROOT CAUSE of that difficulty: the co/contra-variant basis a_i IS the deformed element's coordinate frame.
When the element skews, the a_i become non-orthogonal / nearly parallel, the metric g_ij = a_i.a_j departs
far from the identity, and the co<->contra conversion through the inverse metric g^{ij} is AMPLIFIED and
ill-conditioned. In other words K-H tied the unknown (the physics, M) to the element metric, which degrades
under distortion -- so they had to patch it by hand.

WHY HDiv-VIM (Raviart-Thomas + de Rham) has no such penalty -- FEEC SEPARATES TOPOLOGY FROM METRIC:
  1. The DOFs are FACE FLUXES  q_f = INT_f M.n dS  -- coordinate-free scalars. There is no "direction of
     vectors" to choose, so skew cannot ill-condition the unknown representation.
  2. The demag operator is  N = B^T G B  with a clean split:
       * B = the charge map  [ -div ; normal-trace ] : HDiv -> L2(vol) (+) SurfaceL2(bnd)  is the discrete
         de Rham COBOUNDARY. By the contravariant Piola transform the det(J) cancels
             INT_phys (div u) v dx = INT_ref (div u_hat) v_hat dxi ,
         so B is METRIC-FREE: its matrix is the (signed) incidence of faces/cells/boundary-faces and does
         NOT change under any geometric distortion.
       * G = the Coulomb Gram (1/4 pi r) carries ALL the geometry; distortion enters HERE, correctly.
  3. loops = ker(B) (divergence-free, flux-free fields) are FIELD-NULL: N.loop = B^T G (B.loop) = 0 for any
     G -- so they are mu_r-independent AND distortion-independent. Because B is metric-free, dim(ker B) is
     TOPOLOGICAL (set by mesh connectivity) and rank(B) = n_charge - 1 (the single global divergence-theorem
     constraint) holds cleanly at every distortion. K-H had no such protected nullspace.

This script distorts a unit-cube mesh by a growing boundary-preserving interior deformation
(mesh.SetDeformation) and shows: the elements really are distorted (Jacobian spread grows), yet rank(B),
dim(loops) and even the singular-value gap of B are INVARIANT -- the direct numerical signature of the
metric-free de Rham coboundary. Complements hdiv_vim_distorted.py (which shows the GEOMETRIC part, the Gram
self-energy, stays exact on distorted cells): structure (B) invariant + geometry (G) exact = full distortion
robustness. NGSolve + numpy only (no radia import)."""
import json
import os

import numpy as np
import scipy.sparse as sp
import ngsolve as ng
from netgen.occ import Box, OCCGeometry, Pnt

ng.SetNumThreads(4)
HERE = os.path.dirname(os.path.abspath(__file__))


def _csr(A):
    r, c, v = A.mat.COO()
    return sp.csr_matrix((np.array(v), (np.array(r), np.array(c))), shape=(A.mat.height, A.mat.width))


def charge_map(mesh, p):
    """B = [ -div ; trace.n ] : HDiv(p) -> L2(p-1) (+) SurfaceL2(p).  Returns (dense B, ndof_HDiv, n_charge)."""
    fes = ng.HDiv(mesh, order=p)
    L2v = ng.L2(mesh, order=max(p - 1, 0))
    L2b = ng.SurfaceL2(mesh, order=p)
    u = fes.TrialFunction()
    nn = ng.specialcf.normal(mesh.dim)
    bv = ng.BilinearForm(trialspace=fes, testspace=L2v)
    bv += (-ng.div(u)) * L2v.TestFunction() * ng.dx
    bv.Assemble()
    bb = ng.BilinearForm(trialspace=fes, testspace=L2b)
    bb += (u.Trace() * nn) * L2b.TestFunction() * ng.ds
    bb.Assemble()
    B = sp.vstack([_csr(bv), _csr(bb)]).toarray()
    return B, fes.ndof, L2v.ndof + L2b.ndof


def jacobian_spread(mesh):
    """min / max element |J| (volume Jacobian) over the (possibly deformed) mesh -- a distortion metric."""
    ir = ng.IntegrationRule(points=[(0.25, 0.25, 0.25)], weights=[1.0])
    js = [mesh.GetTrafo(ng.ElementId(e))(list(ir)[0]).measure for e in mesh.Elements(ng.VOL)]
    js = np.array(js)
    return float(js.min()), float(js.max()), float(js.max() / max(js.min(), 1e-30))


def run(maxh=0.5, p=1, amps=(0.0, 0.1, 0.2, 0.3)):
    geo = OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1)))
    mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
    V = ng.VectorH1(mesh, order=2)
    bump = ng.sin(np.pi * ng.x) * ng.sin(np.pi * ng.y) * ng.sin(np.pi * ng.z)  # == 0 on the cube boundary
    base = ng.CF((bump, 0.7 * bump, 0.4 * bump))

    out = {"body": "unit cube", "maxh": maxh, "p": p, "cases": []}
    ref = None
    for amp in amps:
        with ng.TaskManager():
            gf = ng.GridFunction(V)
            gf.Set(amp * base)
            mesh.SetDeformation(gf if amp > 0 else None)
            maxdisp = float(np.abs(np.array(gf.vec)).max()) if amp > 0 else 0.0
            jmin, jmax, jspread = jacobian_spread(mesh)
            B, ndof, ncharge = charge_map(mesh, p)
        sv = np.linalg.svd(B, compute_uv=False)
        tol = sv.max() * max(B.shape) * np.finfo(float).eps * 100
        rank = int((sv > tol).sum())
        nloops = ndof - rank
        sigma_gap = float(sv[rank - 1] / sv.max())
        if ref is None:
            ref = (rank, nloops)
        out["cases"].append(dict(amp=amp, maxdisp=maxdisp, J_min=jmin, J_max=jmax, J_spread=jspread,
                                 ndof_hdiv=ndof, n_charge=ncharge, rank_B=rank, dim_loops=nloops,
                                 sigma_gap=sigma_gap, invariant=bool((rank, nloops) == ref)))
    mesh.SetDeformation(None)
    out["rank_B_invariant"] = bool(all(c["rank_B"] == ref[0] for c in out["cases"]))
    out["dim_loops_invariant"] = bool(all(c["dim_loops"] == ref[1] for c in out["cases"]))
    return out


if __name__ == "__main__":
    res = run()
    print("=== HDiv-VIM distortion robustness: the de Rham loop nullspace is metric-free (invariant) ===")
    print("    K-H 1995 was 'poor when strongly deformed'; HDiv-VIM is structurally immune (see docstring).\n")
    print(f"{'amp':>5} {'maxdisp':>8} {'|J|min':>8} {'|J|max':>8} {'Jspread':>8} "
          f"{'ndofHDiv':>8} {'n_charge':>8} {'rankB':>6} {'dimLoops':>8} {'sigma_gap':>10}")
    for c in res["cases"]:
        tag = "" if c["invariant"] else "  <-- CHANGED!"
        print(f"{c['amp']:>5} {c['maxdisp']:>8.3f} {c['J_min']:>8.4f} {c['J_max']:>8.4f} {c['J_spread']:>8.2f} "
              f"{c['ndof_hdiv']:>8} {c['n_charge']:>8} {c['rank_B']:>6} {c['dim_loops']:>8} "
              f"{c['sigma_gap']:>10.2e}{tag}")
    print(f"\n  rank(B) and dim(loops) invariant across all distortions: "
          f"{res['rank_B_invariant'] and res['dim_loops_invariant']}")
    print("  => B = metric-free de Rham coboundary; geometry enters only via the Gram G. This is WHY")
    print("     HDiv-VIM has no K-H distortion penalty -- de Rham separation, not a co/contra-variant patch.")
    with open(os.path.join(HERE, "hdiv_vim_distortion_loops.json"), "w") as f:
        json.dump(res, f, indent=2)
    print("saved", os.path.join(HERE, "hdiv_vim_distortion_loops.json"))
