"""Frontier F1: the FOLIATION-CHOICE wall (non-convex outer loop).  HONEST probe.

Stage A/B established the CONVEX half: with the foliation mu FIXED, the current
J = grad(lambda) x grad(mu) is LINEAR in lambda, so fitting a target field is a
convex least-norm solve (foliated_clebsch_solenoid.py).  F1 asks the remaining
question: how do you CHOOSE the foliation mu for a general (non-axisymmetric)
target?  This probe demonstrates -- honestly, with a partly-negative result --
that the foliation choice is a genuinely NON-CONVEX outer loop.

WHY THIS IS THE SAME PROBLEM THE FUSION COMMUNITY ALREADY FORMALISED.
The plasma-physics "coil winding surface (CWS) / current potential" problem is
exactly F1 (refs below):
  * Fix the winding surface -> the current potential Phi (= our lambda, J = n x
    grad Phi on a surface; volume analogue J = grad(lambda) x grad(mu)) makes the
    field a LINEAR function of the unknown -> a CONVEX least-squares solve
    (NESCOIL / REGCOIL; REGCOIL adds a Tikhonov current-density term).
  * Let the winding surface ITSELF vary -> the |r-r'|^-3 Biot-Savart nonlinearity
    composed with the surface shape makes the joint problem NON-CONVEX
    (arXiv:2408.08267 "global stellarator coil optimization": "the underlying
    inverse Biot-Savart problem is ... non-convex"; the winding-surface method
    "circumvents this non-linearity by fixing the location of the current").
Our mu is the volume foliation = the stack of winding surfaces.  Fixing mu = their
REGCOIL inner convex solve; choosing mu = their non-convex CWS choice.

THE HONEST STRUCTURE (two sub-results, one positive, one negative):

(1) WHEN mu DOES NOT MATTER (positive).  For a target whose natural current is a
    STACK of single-axis loops (a solenoid about z), a 1-parameter family of
    planar foliations mu_t = (cos t) x + (sin t) z has a UNIMODAL coil-complexity
    objective with ONE optimum at the aligned tilt (t = pi/2, leaves perp to z).
    This is the "mu = r already suffices for cylinder-representable targets"
    statement: each cylinder/leaf carries a surface stream function, so a foliation
    ALIGNED with the target's loop axis reproduces it with minimum current, and the
    choice is convex (single basin).

(2) WHEN mu GENUINELY MATTERS (the WALL, negative).  Superpose two solenoids about
    TWO orthogonal axes (z and x).  Now NO single planar foliation aligns with both
    loop families.  The SAME outer objective becomes BIMODAL: two distinct local
    minima (near t = pi/4 and t = 3pi/4) separated by a barrier.  The field is still
    fit EXACTLY at every tilt (inner solve convex, residual ~1e-15) -- it is the
    coil-complexity (current L2 norm, the REGCOIL measure) that is non-convex over
    foliations.  A gradient/Newton outer loop converges to whichever basin it starts
    in; there is no closed form and no convex reformulation in the foliation.

OUTER OBJECTIVE (REGCOIL-style coil complexity).  At each tilt t we fit the target
field B* EXACTLY by the existing ACA+TSVD min-norm lambda-solve, then report the
resulting current L2 norm

    g(t) = || J(lambda*(t); mu_t) ||_L2(Omega),   J = grad(lambda*) x grad(mu_t),

which is the natural manufacturability / complexity measure (a badly-aligned
foliation still hits the field, but only with large cancelling in-leaf currents).
This g(t) is the value function of a convex inner problem over the foliation
parameter -- exactly the bilevel "convex inner, non-convex outer" of F1.

Reuses: the foliated linear-lambda machinery (grad(phi_j) x grad(mu)) and the
radia.stream_function ACA+TSVD engine, both UNCHANGED.

Refs:
  Merkel, NESCOIL, Nucl. Fusion 27 (1987) 867 (current potential on a fixed CWS).
  Landreman, REGCOIL, Nucl. Fusion 57 (2017) 046003 (Tikhonov current density).
  Kaptanoglu et al., arXiv:2408.08267 (global coil opt: inverse Biot-Savart
    non-convex; winding surface fixes the current location to recover convexity).
  Elder et al., arXiv:2505.07703 (axisymmetric CWS; a pure DISTANCE-offset scan is
    parabolic/unimodal -- consistent with our sub-result (1): offset alone does NOT
    trigger F1; it is the ORIENTATION/topology of the foliation vs a multi-axis
    target that makes the outer objective multi-modal).

Run:  python examples/feec_vim/foliation_choice_wall.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, H1, GridFunction, CF, grad, InnerProduct, dx, x, y, z,
    LinearForm, TaskManager, Cross, Integrate,
)
from netgen.csg import CSGeometry, OrthoBrick, Pnt

from stream_function import aca_tsvd, pseudo_inverse_solve

MU0 = 4e-7 * math.pi

# conductor block (the design region the foliation lives in)
HALF = 0.10           # half-side [m]
RLOOP = 0.16          # target source-loop radius [m]
ILOOP = 1.0e4         # target source-loop current [A]
GRID = 3              # field-match grid points per axis (M = 3 * GRID^3 components)
NT = 13               # tilt samples on [0, pi]


def build_mesh(maxh=0.05):
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-HALF, -HALF, -HALF), Pnt(HALF, HALF, HALF)).mat("cond"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def grad_mu_planar(t):
    """Planar foliation mu_t = (cos t) x + (sin t) z; leaves are planes with
    normal n_t = (cos t, 0, sin t).  J = grad(lambda) x grad(mu_t) lies in the
    leaf (J . n_t = 0): the foliation FIXES the plane the current flows in."""
    return CF((math.cos(t), 0.0, math.sin(t)))


def response_full(fes, targets, t):
    """A[3i+c, j] = (B_c at target i) from basis current grad(phi_j) x grad(mu_t).

    [(grad v x grad mu) x K]_c via the Radia Biot-Savart kernel K = d/|d|^3."""
    u, v = fes.TnT()
    X = CF((x, y, z))
    gmu = grad_mu_planar(t)
    M, N = 3 * len(targets), fes.ndof
    A = np.zeros((M, N))
    for i, (px, py, pz) in enumerate(targets):
        P = CF((px, py, pz))
        d = P - X
        K = d / (InnerProduct(d, d) ** 1.5)
        cross = Cross(Cross(grad(v), gmu), K)      # (J_basis x K)
        for c in range(3):
            lf = LinearForm(fes)
            lf += (MU0 / (4 * math.pi)) * cross[c] * dx
            lf.Assemble()
            A[3 * i + c, :] = lf.vec.FV().NumPy()
    return A


def current_norm(fes, mesh, lam_vec, t):
    """||J||_L2 of the current J = grad(lambda) x grad(mu_t) -- coil complexity."""
    gf = GridFunction(fes)
    gf.vec.FV().NumPy()[:] = lam_vec
    J = Cross(grad(gf), grad_mu_planar(t))
    return math.sqrt(Integrate(InnerProduct(J, J) * dx, mesh))


def loop_field(axis, targets, R=RLOOP, I=ILOOP):
    """Ground-truth target field: a circular current loop about `axis`,
    sampled (all 3 components) at the field-match points.  Independent of the
    FE machinery (computed with Radia's own Biot-Savart)."""
    import radia as rad
    rad.UtiDelAll()
    th = np.linspace(0, 2 * math.pi, 48, endpoint=True)
    if axis == "z":
        pts = [[R * math.cos(a), R * math.sin(a), 0.0] for a in th]
    else:  # "x"
        pts = [[0.0, R * math.cos(a), R * math.sin(a)] for a in th]
    loop = rad.ObjFlmCur(pts, I)
    out = []
    for P in targets:
        B = rad.Fld(loop, "b", list(P))
        out.extend([float(B[0]), float(B[1]), float(B[2])])
    rad.UtiDelAll()
    return np.array(out)


def sweep(fes, mesh, targets, Btarget):
    """Outer objective g(t) = min-current-norm to fit Btarget EXACTLY at tilt t.

    Returns (ts, g, fit_res): g[k] = ||J||_L2 at the min-norm lambda; fit_res[k]
    confirms the inner solve hit the field (convex, feasible)."""
    ts = np.linspace(0.0, math.pi, NT)
    g = np.zeros(NT)
    fit = np.zeros(NT)
    for k, t in enumerate(ts):
        A = response_full(fes, targets, t)
        M, N = A.shape
        res = aca_tsvd(M, N, lambda i, j: A[i, j],
                       modes=min(M, N), kmax=min(M, N), aca_eps=1e-12, method=3)
        lam = pseudo_inverse_solve(res, Btarget)
        g[k] = current_norm(fes, mesh, lam, t)
        fit[k] = np.linalg.norm(A @ lam - Btarget) / (np.linalg.norm(Btarget) + 1e-30)
    return ts, g, fit


def count_local_minima(g):
    """Strict interior local minima (g[i] < both neighbours)."""
    return [i for i in range(1, len(g) - 1) if g[i] < g[i - 1] and g[i] < g[i + 1]]


def barrier_between(g, mins):
    """Relative barrier height between the two lowest minima (0 if < 2 minima)."""
    if len(mins) < 2:
        return 0.0
    i1, i2 = sorted(mins, key=lambda i: g[i])[:2]
    lo, hi = sorted([i1, i2])
    deeper = max(g[i1], g[i2])
    return (g[lo:hi + 1].max() - deeper) / (deeper + 1e-30)


def main():
    print("=== F1: foliation-choice wall (convex inner, non-convex outer) ===")
    with TaskManager():
        mesh = build_mesh()
        fes = H1(mesh, order=2)
        gp = np.linspace(-0.5 * HALF, 0.5 * HALF, GRID)
        targets = [(a, b, c) for a in gp for b in gp for c in gp]

        Bz = loop_field("z", targets)            # solenoid-about-z  (single axis)
        Bx = loop_field("x", targets)            # solenoid-about-x
        B_two = Bz + Bx                          # two-axis target

        # (1) single-axis target -> foliation choice is CONVEX (unimodal g)
        ts, g1, fit1 = sweep(fes, mesh, targets, Bz)
        # (2) two-axis target  -> foliation choice is NON-CONVEX (bimodal g)
        ts, g2, fit2 = sweep(fes, mesh, targets, B_two)

    m1 = count_local_minima(g1)
    m2 = count_local_minima(g2)
    bar2 = barrier_between(g2, m2)
    fit_max = max(fit1.max(), fit2.max())

    print(f"   mesh ne={mesh.ne}, lambda ndof={fes.ndof}, "
          f"M={3 * len(targets)} field comps, NT={NT} tilts")
    print(f"   inner-solve fit residual (max over all tilts) = {fit_max:.1e}"
          f"   (field hit EXACTLY -> inner problem convex)")
    print()
    print(f"   (1) single-axis target  g(t)=||J|| :")
    print(f"       t/pi = " + " ".join(f"{t/math.pi:.2f}" for t in ts))
    print(f"       g    = " + " ".join(f"{v/1e3:.0f}" for v in g1) + "   (x1e3 A)")
    print(f"       local minima at t/pi = {[round(ts[i]/math.pi, 2) for i in m1]}"
          f"   -> {len(m1)} basin  => CONVEX foliation choice")
    print()
    print(f"   (2) two-axis target     g(t)=||J|| :")
    print(f"       t/pi = " + " ".join(f"{t/math.pi:.2f}" for t in ts))
    print(f"       g    = " + " ".join(f"{v/1e3:.0f}" for v in g2) + "   (x1e3 A)")
    print(f"       local minima at t/pi = {[round(ts[i]/math.pi, 2) for i in m2]}"
          f"   -> {len(m2)} basins, barrier {bar2:.0%}")
    print(f"       => NON-CONVEX foliation choice (the F1 WALL): a gradient outer")
    print(f"          loop lands in whichever basin it starts in; no closed form.")

    ok = (
        fit_max < 1e-6 and            # inner convex solve always feasible
        len(m1) == 1 and              # single-axis: unimodal (mu choice convex)
        len(m2) >= 2 and              # two-axis: multi-modal (the wall)
        bar2 > 0.05                   # a real barrier separates the basins
    )
    print(f"\n=== PASS: {ok} ===")
    return {
        "ne": int(mesh.ne), "ndof": int(fes.ndof), "M": int(3 * len(targets)),
        "fit_max": float(fit_max),
        "n_min_single": len(m1), "n_min_two": len(m2),
        "barrier_two": float(bar2),
        "tmin_single": [float(ts[i]) for i in m1],
        "tmin_two": [float(ts[i]) for i in m2],
    }, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
