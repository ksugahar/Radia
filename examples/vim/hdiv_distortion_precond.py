"""
hdiv_distortion_precond.py -- the LAST open piece for "accurate + loop-free + solvable-when-distorted"
HDiv-VIM demag: does a METRIC preconditioner bound the star-block conditioning (and the MINRES iters)
as the mesh DISTORTION grows?  Everything else is already demonstrated in the sibling prototypes:

  - accuracy        : hdiv_demag_quad_self.py  -> geometry-exact subdivided Gram (vol-err 1e-11 regular,
                      2.6e-3 at 18% distortion; -> exact with the C++ Wilton analytic face integral).
  - loop-free       : loop_res ~4e-16 on regular AND distorted (structural: loops = ker(charge map)).
  - mu_r-independent: MINRES iters 8/8/8 (regular), 80/76/69 (distorted) -- FLAT in mu_r (loops sit at
                      EXACTLY 1/chi, never pollute).  No loop-star, no H-ILU needed for high mu_r.

The ONLY thing that grows is iters-vs-DISTORTION with PLAIN POINT-Jacobi (8 -> ~80 over 0..24%): a
mesh-quality conditioning effect, NOT a mu_r blow-up.  yano_star_scalability.py showed the SAME demag
operator, on the STAR space, with the cell-graph-Laplacian basis-Gram metric M = B_star^T B_star,
conditions to cond ~9..14 (vs raw ~300..18000) at FIXED distortion across sizes.  Nobody has yet swept
INCREASING distortion with that metric.  This script does exactly that:

  for distortion 0..0.24 (mu_r=1e4):
    A      = (1/chi) I - N            (full HDiv-dof system; N = B^T G B, geometry-exact G)
    A_star = (1/chi) I - Nsb          (SVD-orthonormal star; loops removed)
    report cond + MINRES iters for:
      (P0) point |diag| Jacobi on the FULL system          -- the baseline that grows 8->80
      (P1) UNpreconditioned star solve                     -- removing loops alone enough?
      (P2) cell-graph-Laplacian metric on the star         -- M_lap = (B Sstar)^T (B Sstar)
      (P3) HDiv-mass metric on the star                    -- M_m  = Sstar^T mass_hdiv Sstar

MEASURED RESULT (2026-06-21) -- the metric-preconditioner hypothesis is REFUTED, and the real
limit is identified:
  - NO preconditioner bounds the iters across distortion.  Over 0..24%: P0 point-Jacobi 8->82
    (10x), star-raw 36->85, star+H(div)ip 46->126, star+HDivmass 31->69, star+TRUE-cell-graph-
    Laplacian 53->169 (all GROW; the Laplacian metric -- which bounded the SIZE axis in
    yano_star_scalability -- does NOT bound the DISTORTION axis).
  - WHY: the diagnostic gap = mu_min_star * chi collapses 212 -> 0.84 at 24% distortion.  A spurious
    near-zero-demag STAR mode (a bad-element discretisation artifact) descends toward the loop
    position 1/chi, so A_star = (1/chi)I - N_star becomes near-SINGULAR on that mode.  That is an
    OPERATOR near-singularity (mesh-quality), which NO SPD metric preconditioner can remove -- the
    same wall any FEM hits with badly-shaped elements.  It is also a mu_r-position resonance (worst
    when 1/chi lands on the spurious small eigenvalue), not a monotone mu_r blow-up.
  - WHAT HOLDS: at PRACTICAL distortion (<=18%) the gap stays 60..210 (safe) and the system solves
    in ~70..85 MINRES iters, mu_r-INDEPENDENT, loop_res ~4e-16.  So "accurate(geom-exact Gram) +
    loop-free(structural) + solvable-when-distorted" is ALREADY delivered in the usable regime; the
    distortion ceiling is a mesh-quality property, not a solver/preconditioner deficiency.
  => Production recipe (the real gap = it lives only in these Python prototypes): port the exact
     Wilton analytic charged-face Gram (accuracy; the same charged-face field kernel the yano C++
     port exercised) + a MINRES (symmetric-indefinite) solve into the C++ HDiv-VIM, which today
     still uses the crude centroid-monopole Gram.  A distortion-optimal preconditioner is NOT the
     missing piece; mesh quality is the ceiling.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from hdiv_demag_quad_self import build  # geometry-exact Gram builder (imports ngsolve; fast)
from scipy.sparse.linalg import minres, LinearOperator


def cond_precond(A, M):
    """cond of the symmetrically-preconditioned operator M^{-1/2} A M^{-1/2} (M SPD)."""
    w, V = np.linalg.eigh(0.5 * (M + M.T))
    w = np.clip(w, 1e-30, None)
    Mihalf = (V / np.sqrt(w)) @ V.T
    Ap = Mihalf @ A @ Mihalf
    return float(np.linalg.cond(Ap)), Mihalf


def minres_iters(A, M_inv_op):
    n = A.shape[0]
    it = {"n": 0}
    minres(A, np.ones(n), M=M_inv_op, rtol=1e-8, maxiter=6000,
           callback=lambda xk: it.__setitem__("n", it["n"] + 1))
    return it["n"]


def spd_inv_op(M):
    w, V = np.linalg.eigh(0.5 * (M + M.T))
    w = np.clip(w, 1e-30, None)
    Minv = (V / w) @ V.T
    return LinearOperator(M.shape, matvec=lambda x: Minv @ x)


def run():
    MUR = 1e4
    chi = MUR - 1.0
    sweep = []
    print("=== HDiv-VIM demag: does a star METRIC preconditioner bound conditioning vs distortion? ===")
    print(f"    (mu_r={MUR:.0e}, geometry-exact Gram, 3x3x3 hexes, ndof=108, star=80, loops=28)\n")
    hdr = (f"{'distort':>7} | {'loop_res':>9} | {'P0 pointJac':>11} | {'P1 star raw':>15} | "
           f"{'P2 H(div)ip':>15} | {'P3 HDivmass':>15} | {'P4 Laplacian':>15} | {'gap':>10}")
    print(hdr)
    print("-" * len(hdr))
    for da in (0.0, 0.06, 0.12, 0.18, 0.24):
        d = build(da, 6, 8)
        ndof, rankQ = d["ndof"], d["rankQ"]
        N, Nsb, Sstar = d["N"], d["Nsb"], d["Sstar"]
        # full system + star system
        A = (1.0 / chi) * np.eye(ndof) - N
        A_star = (1.0 / chi) * np.eye(rankQ) - Nsb
        # SPD metrics on the star space (loops already removed by the SVD projection):
        #   P3 = HDiv L2 mass; P2 = H(div) inner product (mass + div-div), the AMS-type metric;
        #   P4 = the TRUE cell-graph Laplacian B_star^T B_star (the metric that bounded the SIZE
        #        axis in yano_star_scalability.py -- B = the charge map, now exposed by build()).
        M_mass_star = Sstar.T @ d["mass_hdiv"] @ Sstar
        Hdiv_ip = d["mass_hdiv"] + d["divdiv_hdiv"]
        M_hdivip_star = Sstar.T @ Hdiv_ip @ Sstar
        Bstar = d["B"] @ Sstar
        M_lap_star = Bstar.T @ Bstar

        # P0: point |diag| Jacobi on the FULL system (the baseline that grows)
        dgA = np.abs(np.diag(A)).copy(); dgA[dgA < 1e-30] = 1.0
        it_p0 = minres_iters(A, LinearOperator(A.shape, matvec=lambda x, dgA=dgA: x / dgA))
        cond_p0 = float(np.linalg.cond(A))

        # P1: star, unpreconditioned (identity)
        it_p1 = minres_iters(A_star, LinearOperator(A_star.shape, matvec=lambda x: x))
        cond_p1 = float(np.linalg.cond(A_star))

        # P2: star + H(div) inner-product metric (mass + div-div) -- the AMS-type metric
        cond_p2, _ = cond_precond(A_star, M_hdivip_star)
        it_p2 = minres_iters(A_star, spd_inv_op(M_hdivip_star))

        # P3: star + HDiv mass metric
        cond_p3, _ = cond_precond(A_star, M_mass_star)
        it_p3 = minres_iters(A_star, spd_inv_op(M_mass_star))

        # P4: star + TRUE cell-graph Laplacian metric (B_star^T B_star)
        cond_p4, _ = cond_precond(A_star, M_lap_star)
        it_p4 = minres_iters(A_star, spd_inv_op(M_lap_star))

        # near-null diagnostic: how close is the smallest star demag mode to the loop position 1/chi?
        gap = float(d["mu"][0] * chi)  # mu_min_star / (1/chi); ~1 => a star mode resonates with loops

        print(f"{da:>7.2f} | {d['loop_res']:>9.1e} | {it_p0:>4d} (c{cond_p0:>7.0f}) | "
              f"{it_p1:>4d} ({cond_p1:>8.1f}) | {it_p2:>4d} ({cond_p2:>8.1f}) | "
              f"{it_p3:>4d} ({cond_p3:>8.1f}) | {it_p4:>4d} ({cond_p4:>8.1f}) | gap={gap:>6.2f}")
        sweep.append({
            "distort": da, "loop_res": float(d["loop_res"]), "mu_min_star": float(d["mu"][0]),
            "star_loop_gap": gap,
            "P0_pointjac": {"iters": int(it_p0), "cond_full": cond_p0},
            "P1_star_raw": {"iters": int(it_p1), "cond": cond_p1},
            "P2_star_hdiv_ip": {"iters": int(it_p2), "cond": cond_p2},
            "P3_star_hdiv_mass": {"iters": int(it_p3), "cond": cond_p3},
            "P4_star_cellgraph_laplacian": {"iters": int(it_p4), "cond": cond_p4},
        })
    # data-driven verdict (no pre-judging): does ANY preconditioner keep iters within 2x of its
    # distort=0 value across the sweep?
    print("\nVerdict (data-driven):")
    keys = [("P0_pointjac", "point-Jacobi(full)"), ("P1_star_raw", "star raw"),
            ("P2_star_hdiv_ip", "star+H(div)ip"), ("P3_star_hdiv_mass", "star+HDivmass"),
            ("P4_star_cellgraph_laplacian", "star+cell-graph Laplacian")]
    for k, name in keys:
        it0 = sweep[0][k]["iters"]; itN = sweep[-1][k]["iters"]
        ratio = itN / max(it0, 1)
        flag = "BOUNDED" if ratio <= 2.0 else "GROWS"
        print(f"  {name:<28}: iters {it0:>4d} -> {itN:>4d}  ({ratio:.1f}x)  {flag}")
    gaps = [s["star_loop_gap"] for s in sweep]
    print(f"  star-loop gap (mu_min_star*chi) over sweep: {[round(g,2) for g in gaps]}")
    print("  (gap -> ~0 means a star demag mode collapses toward the loop position 1/chi: an")
    print("   OPERATOR near-singularity that NO SPD metric preconditioner can remove.)")
    out = {"mu_r": MUR, "grid": "3x3x3 hexes", "ndof": 108, "sweep": sweep}
    with open(os.path.join(HERE, "hdiv_distortion_precond.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nsaved", os.path.join(HERE, "hdiv_distortion_precond.json"))


if __name__ == "__main__":
    run()
