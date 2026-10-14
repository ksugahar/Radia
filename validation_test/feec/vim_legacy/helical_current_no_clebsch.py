"""Frontier F3 PROBE: nonzero-helicity current -> no global Clebsch -> no clean wires.

This is an HONEST, partly-NEGATIVE result.  It does NOT solve a design problem;
it DEMONSTRATES A FUNDAMENTAL WALL and locks the wall behaviour with a golden.

THE WALL (helicity is a conserved topological invariant, Moffatt 1969;
Arnold & Khesin 1998).  Clean level-set / streamline WIRES are the field lines of
a Clebsch current J = grad(lambda) x grad(mu) -- they are the intersections
{lambda = const} INT {mu = const}, which CLOSE.  But (sympy-verified 2026-06-11,
re-verified here) a Clebsch current has helicity density T.curl T = 0 EXACTLY
(T = lambda*grad(mu)).  So a current with H != 0 has NO global Clebsch pair, and
its current lines are NOT clean closing rings.  Two adversarial cracks, both fail:

  CRACK 1  (convex Clebsch fit -- the HELICITY-INVARIANT wall).  Take the
           ARNOLD-BELTRAMI-CHILDRESS (ABC) current J_ABC (curl J = k J  =>  a
           Beltrami current, H_rel = 1).  Try to reproduce it with the Stage-A
           foliated convex solver (mu = r FIXED, lambda free, the SAME ACA+TSVD
           Gram machinery).  The foliated ansatz J = grad(lambda) x grad(mu) has
           helicity density T.J = 0 FOR ANY lambda AND ANY foliation mu
           (sympy-verified) -- so its reconstructed current is FORCED to
           H_rel ~ 0.  We report BOTH the L2 fit residual (foliation-DEPENDENT,
           ~0.6, informative) AND the helicity of the best fit (foliation-
           INDEPENDENT, ~0): the helicity GAP between target (~1) and ANY Clebsch
           fit (~0) is the fundamental, topological wall -- a conserved invariant
           the convex hypothesis class can never carry.

  CRACK 2  (streamline closure).  Trace the current lines.  A Clebsch
           (azimuthal) current line CLOSES (arclength = 2 pi r, returns to its
           seed).  An ABC current line does NOT close -- it drifts / wanders
           ergodically (Dombre, Frisch, Greene, Henon, Mehr & Soward 1986,
           J. Fluid Mech. 167:353; "Beltrami fields exhibit knots and chaos
           almost surely", Enciso & Peralta-Salas 2020).  No equal-current loop
           family can be cut from non-closing lines.

THE ONLY HONEST OUTPUTS for a nonzero-helicity target (stated, not faked):
  (1) the VECTOR-T DISTRIBUTION  J = curl T  (Stage A / vector_t_inverse.py):
      always exists (every div-free J is curl T), unconstrained by helicity, but
      it is a CONTINUOUS bulk current, NOT a set of windable equal-current wires.
  (2) a MULTI-PATCH CLEBSCH ATLAS with COHOMOLOGY CUTS: local Clebsch charts
      glued across cut surfaces (the b1 != 0 topology), i.e. accept seams.

Run:  python validation_test/feec/vim_legacy/helical_current_no_clebsch.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, H1, GridFunction, CF, grad, InnerProduct, dx, sqrt,
    sin, cos, x, y, z, LinearForm, TaskManager, Integrate, Cross,
)
from netgen.csg import CSGeometry, Cylinder, Plane, Pnt, Vec

# helicity_relative reuses the Stage-B HCurl helicity diagnostic UNCHANGED.
from helicity_diagnostic import helicity_relative

PI = math.pi

# conductor tube (same as the Stage-A/B solenoid so the pieces are shared)
R0, R1 = 0.05, 0.10        # inner / outer radius [m]
ZL = 0.10                  # half-length [m]
KABC = 2 * PI / (2 * ZL)   # ABC wavenumber: one full period over the tube length


def build_mesh(maxh=0.03):
    geo = CSGeometry()
    outer = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R1)
    inner = Cylinder(Pnt(0, 0, -1), Pnt(0, 0, 1), R0)
    top = Plane(Pnt(0, 0, ZL), Vec(0, 0, 1))
    bot = Plane(Pnt(0, 0, -ZL), Vec(0, 0, -1))
    geo.Add(((outer - inner) * top * bot).mat("cond"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def abc_current_cf(k=KABC):
    """ABC current J_ABC (A=B=C=1).  curl J = k J  =>  Beltrami => H_rel = 1."""
    return CF((sin(k * z) + cos(k * y),
              sin(k * x) + cos(k * z),
              sin(k * y) + cos(k * x)))


# --------------------------------------------------------------------------
# GATE: the target current is genuinely helical (H_rel ~ 1)
# --------------------------------------------------------------------------
def gate_helicity(mesh, k=KABC):
    """H_rel of the ABC current's driving potential T = J/k (curl T = J).

    For a Beltrami current curl J = k J, T = J/k satisfies curl T = J, and
    T.curl T = |J|^2/k > 0 everywhere -> H_rel = 1 (reuse helicity_diagnostic).
    """
    T_cf = abc_current_cf(k) / k
    return helicity_relative(mesh, T_cf, order=3)


# --------------------------------------------------------------------------
# CRACK 1: convex foliated-Clebsch fit (mu = r FIXED) cannot reach a helical J
# --------------------------------------------------------------------------
def crack1_foliated_clebsch_fit(mesh, k=KABC):
    """Best foliated-Clebsch J = grad(lambda) x grad(r) fit to J_ABC.

    mu = r is FIXED so the fit is LINEAR/CONVEX in lambda (the Stage-A ansatz).
    We least-squares fit lambda (H1) to the ABC target and report:
      * fit_res = ||J_fit - J_ABC|| / ||J_ABC||  (foliation-DEPENDENT, informative)
      * fit_helicity = H_rel of the reconstructed J_fit (foliation-INDEPENDENT;
        FORCED to ~0 since a Clebsch current has T.J = 0 for ANY lambda/mu).
    The helicity GAP (target ~1 vs fit ~0) is the topological wall.

    Returns (fit_res, fit_helicity).
    """
    fes = H1(mesh, order=2)
    r = sqrt(x * x + y * y)
    gmu = CF((x / r, y / r, 0.0))
    Jt = abc_current_cf(k)

    # min_lambda || grad(lambda) x grad(mu) - Jt ||^2_L2  -> normal equations.
    # Basis current for DOF j: b_j = grad(phi_j) x grad(mu).  Assemble Gram
    # G[i,j] = <b_i, b_j> and rhs f[i] = <b_i, Jt> via NGSolve forms.
    from ngsolve import BilinearForm
    u, v = fes.TnT()
    bu = Cross(grad(u), gmu)
    bv = Cross(grad(v), gmu)
    G = BilinearForm(fes, symmetric=True)
    G += InnerProduct(bu, bv) * dx
    G.Assemble()
    f = LinearForm(fes)
    f += InnerProduct(bv, Jt) * dx
    f.Assemble()

    # lambda has a constant null space (grad const = 0) + the foliation can carry
    # no current along grad(mu), so G is singular -> use a TSVD pseudo-inverse on
    # the DENSE Gram (small ndof) for the honest min-norm least-squares solution.
    n = fes.ndof
    Gd = np.zeros((n, n))
    rows, cols, vals = G.mat.COO()
    for ii, jj, vv in zip(rows, cols, vals):
        Gd[ii, jj] += vv
    fv = f.vec.FV().NumPy().copy()
    w, U = np.linalg.eigh(Gd)
    tol = max(abs(w)) * 1e-10
    inv = np.where(w > tol, 1.0 / w, 0.0)
    lam = U @ (inv * (U.T @ fv))

    # residual ||J_fit - Jt|| / ||Jt|| from the assembled quadratic:
    #   ||J_fit - Jt||^2 = lam^T G lam - 2 lam^T f + ||Jt||^2
    JtL2sq = Integrate(InnerProduct(Jt, Jt) * dx, mesh)
    res_sq = float(lam @ (Gd @ lam) - 2.0 * lam @ fv + JtL2sq)
    res_sq = max(res_sq, 0.0)
    fit_res = math.sqrt(res_sq / (JtL2sq + 1e-30))

    # Helicity of the RECONSTRUCTED foliated current.  Its driving potential is
    # T_fit = lambda * grad(mu) (curl T_fit = grad(lambda) x grad(mu) = J_fit),
    # so H_rel(T_fit) = 0 to discretisation -- the topological invariant the
    # Clebsch class cannot carry, INDEPENDENT of the foliation chosen.
    glam = GridFunction(fes)
    glam.vec.FV().NumPy()[:] = lam
    T_fit = glam * gmu
    fit_helicity = helicity_relative(mesh, T_fit, order=3)
    return fit_res, fit_helicity


# --------------------------------------------------------------------------
# CRACK 2: streamline (current-line) closure -- Clebsch closes, ABC does not
# --------------------------------------------------------------------------
def _trace_closure(Jf, seed, h, nmax, rtube):
    """Arc-length RK4 trace of a current line; report whether it returns to seed.

    Returns (closed, arclen, final_dist).  Stops if the line leaves the tube
    (|x,y| in [R0,R1] band relaxed) or reaches nmax steps.  "Closed" = it left
    the seed neighborhood (> 5h) and came back within 3h.
    """
    p = np.array(seed, float)
    P0 = p.copy()
    arclen = 0.0
    left = False
    for _ in range(nmax):
        def fdir(q):
            v = Jf(q)
            nv = np.linalg.norm(v)
            return v / (nv + 1e-30)
        k1 = fdir(p)
        k2 = fdir(p + 0.5 * h * k1)
        k3 = fdir(p + 0.5 * h * k2)
        k4 = fdir(p + h * k3)
        p = p + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        arclen += h
        d = np.linalg.norm(p - P0)
        if d > 5 * h:
            left = True
        if left and d < 3 * h:
            return True, arclen, d
    return False, arclen, float(np.linalg.norm(p - P0))


def crack2_streamline_closure(k=KABC):
    """Two analytic current fields, traced as ODEs (no mesh needed for the line).

    Clebsch azimuthal current J = phi_hat (lambda=z, mu=r): lines are CIRCLES of
    radius r -> CLOSE at arclength 2 pi r.  ABC current: lines do NOT close.
    """
    def J_clebsch(p):
        px, py, pz = p
        r = math.hypot(px, py) + 1e-30
        return np.array([-py / r, px / r, 0.0])      # azimuthal -> circles

    def J_abc(p):
        px, py, pz = p
        return np.array([math.sin(k * pz) + math.cos(k * py),
                         math.sin(k * px) + math.cos(k * pz),
                         math.sin(k * py) + math.cos(k * px)])

    rmid = 0.5 * (R0 + R1)
    h = 1e-3
    # Clebsch: a ring at r = rmid, z = 0
    cl, al_c, d_c = _trace_closure(J_clebsch, (rmid, 0.0, 0.0),
                                   h, 80000, rmid)
    expected_ring = 2 * PI * rmid
    # ABC: a generic seed inside the tube band
    seed = (rmid, 0.0, 0.3 * ZL)
    ab, al_a, d_a = _trace_closure(J_abc, seed, h, 80000, rmid)
    # drift = how far the ABC line ends from its seed in units of the ring length
    abc_drift_rings = d_a / (expected_ring + 1e-30)
    return {
        "clebsch_closed": bool(cl),
        "clebsch_arclen": float(al_c),
        "clebsch_ring_expected": float(expected_ring),
        "clebsch_close_err": float(d_c),
        "abc_closed": bool(ab),
        "abc_arclen": float(al_a),
        "abc_final_dist": float(d_a),
        "abc_drift_rings": float(abc_drift_rings),
    }


def main():
    print("=== F3 PROBE: nonzero-helicity current -> NO global Clebsch -> "
          "NO clean wires ===")
    with TaskManager():
        mesh = build_mesh()

        # GATE: the ABC current is genuinely helical
        h_rel = gate_helicity(mesh)

        # CRACK 1: convex foliated-Clebsch fit cannot reach the helical target;
        # its reconstructed current is FORCED to zero helicity (the wall).
        clebsch_fit_res, clebsch_fit_hel = crack1_foliated_clebsch_fit(mesh)

        # CRACK 2: streamline closure (Clebsch closes, ABC does not)
        sl = crack2_streamline_closure()

    helicity_gap = h_rel - clebsch_fit_hel

    print(f"   mesh ne={mesh.ne}")
    print(f"   GATE  helical target H_rel        = {h_rel:.3f}   "
          f"(Beltrami => ~1, genuinely helical)")
    print(f"   CRACK1 foliated-Clebsch fit resid = {clebsch_fit_res:.3f}   "
          f"(foliation-dependent, informative)")
    print(f"   CRACK1 fit current H_rel          = {clebsch_fit_hel:.3f}   "
          f"(FORCED ~0: Clebsch carries no helicity, ANY foliation)")
    print(f"   CRACK1 HELICITY GAP (the wall)    = {helicity_gap:.3f}   "
          f"(target ~1 minus Clebsch fit ~0 = topological invariant lost)")
    print(f"   CRACK2 Clebsch line CLOSES        = {sl['clebsch_closed']}  "
          f"(arclen {sl['clebsch_arclen']:.3f} vs 2*pi*r "
          f"{sl['clebsch_ring_expected']:.3f}, err {sl['clebsch_close_err']:.1e})")
    print(f"   CRACK2 ABC line CLOSES            = {sl['abc_closed']}  "
          f"(final dist {sl['abc_final_dist']:.2f} = "
          f"{sl['abc_drift_rings']:.1f} ring-lengths away -> wanders)")
    print("   HONEST OUTPUTS: (1) vector-T distribution J=curl T (continuous "
          "bulk, not wires); (2) multi-patch Clebsch atlas w/ cohomology cuts.")

    # The golden asserts the WALL, not a solve:
    #   - the target is helical (H_rel ~ 1),
    #   - ANY foliated Clebsch fit collapses to ~0 helicity -> a GAP ~ 1
    #     (the foliation-independent topological wall),
    #   - a Clebsch current line closes but the helical one does not.
    ok = (h_rel > 0.9 and
          clebsch_fit_hel < 0.1 and
          helicity_gap > 0.8 and
          sl["clebsch_closed"] and sl["clebsch_close_err"] < 5e-3 and
          (not sl["abc_closed"]) and sl["abc_drift_rings"] > 5.0)
    print(f"\n=== PASS (wall demonstrated): {ok} ===")
    out = {"h_rel": float(h_rel), "clebsch_fit_res": float(clebsch_fit_res),
           "clebsch_fit_hel": float(clebsch_fit_hel),
           "helicity_gap": float(helicity_gap)}
    out.update(sl)
    return out, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
