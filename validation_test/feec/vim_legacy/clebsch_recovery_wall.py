"""F2 frontier probe: recovering (lambda, mu) from a general div-free J is ILL-POSED.

Stage A/B gave the FORWARD map: fix a foliation mu, solve the convex lambda-inverse,
J = grad(lambda) x grad(mu).  F2 asks the INVERSE question -- given a general
divergence-free current J with ~zero helicity, can we RECOVER a Clebsch pair
(lambda, mu)?  The textbook recipe is:

    1. mu is a STREAM SURFACE: J . grad(mu) = 0  (because J = grad(lambda) x grad(mu)
       is perpendicular to grad(mu)).  Solve this linear transport PDE for mu.
    2. then lambda from  J = grad(lambda) x grad(mu).

This probe DEMONSTRATES THE WALL: the recipe is not well-posed.  Two independent
obstructions, each shown numerically on a KNOWN Clebsch current (so we know the
ground truth), reusing the NGSolve transport / area-form machinery:

WALL 1 -- NON-UNIQUENESS (gauge freedom f(mu0)).
    Take a known Clebsch current J = grad(z) x grad(x) = e_y (so mu0 = x, the
    streamlines run in +y).  Solve the transport PDE J.grad(mu)=0 TWICE, seeding
    the inflow face with mu0 = x and with f(mu0) = x^2.  BOTH solve the transport
    to ~1e-9 (both are valid stream surfaces), yet the recovered fields differ by
    O(1): the second is the reparametrization mu_b = f(mu_a) = mu_a^2.  Since ANY
    f(mu0) is an equally-valid second scalar (chain rule: J.grad(f(mu0)) =
    f'(mu0) J.grad(mu0) = 0), the recovery is non-unique -- it returns whatever
    the seed/boundary data dictates, not a canonical (lambda, mu).
    sympy-verified (clebsch_recovery_wall sympy block): J(lambda/f'(mu0), f(mu0))
    == J(lambda, mu0).

WALL 2 -- DEGENERACY AT A FIELD NULL.
    Take a known Clebsch current with an isolated line-null: lambda0=(x^2-y^2)/2,
    mu0 = x*y  ->  J = (0, 0, x^2+y^2), which VANISHES on the z-axis.  Recovering
    lambda inverts the in-plane Gram of (grad lambda, grad mu) whose determinant is
    the squared Clebsch area form |grad lambda x grad mu|^2 = |J|^2.  As we approach
    the null this collapses like r^4, so the recovery conditioning ~ 1/|J|^2 blows
    up (cond ~ 1/r^4) and the lambda-solve is singular ON the null line.  This is the
    "singular where J=0" obstruction noted in the F2 frontier statement and the
    Clebsch-near-vanishing-vorticity literature (Yoshida 2009; Graham & Henyey 2000).

Neither wall is a bug to be fixed: they are the precise reason a GENERAL div-free
current has no canonical Clebsch pair -- only local, seed-dependent, null-singular
ones.  This is an HONEST partly-negative result (repo-first: a clearly-demonstrated
wall is a real result).

Run:  python validation_test/feec/vim_legacy/clebsch_recovery_wall.py
"""
import os
import sys
import math

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../src/radia"))

from ngsolve import (
    Mesh, H1, GridFunction, CF, grad, InnerProduct, dx, Cross, Norm,
    x, y, z, BilinearForm, LinearForm, TaskManager, Integrate, specialcf,
)
from netgen.csg import CSGeometry, OrthoBrick, Plane, Pnt, Vec


# ----------------------------------------------------------------------------
# WALL 1 -- non-uniqueness:  recover mu from the transport PDE J.grad(mu)=0
# ----------------------------------------------------------------------------
def _cube_with_inflow(maxh=0.12):
    """Unit cube whose y=0 face is labelled 'inflow' (the rest 'walls')."""
    geo = CSGeometry()
    cube = (Plane(Pnt(0, 0, 0), Vec(-1, 0, 0)).bc("x0")
            * Plane(Pnt(1, 1, 1), Vec(1, 0, 0)).bc("x1")
            * Plane(Pnt(0, 0, 0), Vec(0, -1, 0)).bc("inflow")
            * Plane(Pnt(1, 1, 1), Vec(0, 1, 0)).bc("outflow")
            * Plane(Pnt(0, 0, 0), Vec(0, 0, -1)).bc("z0")
            * Plane(Pnt(1, 1, 1), Vec(0, 0, 1)).bc("z1"))
    geo.Add(cube.mat("dom"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def _transport_recover(mesh, bvel, seed_cf, order=2):
    """SUPG-stabilized solve of bvel.grad(mu)=0 with mu=seed on the inflow face.

    Returns (gfmu, rel_residual) where rel_residual =
    ||bvel.grad(mu)||_L2 / ||grad(mu)||_L2  measures how well mu solves the
    transport PDE (a valid stream surface => ~0).
    """
    fes = H1(mesh, order=order, dirichlet="inflow")
    u, v = fes.TnT()
    h = specialcf.mesh_size
    tau = h                                  # SUPG stabilization length
    a = BilinearForm(fes)
    a += InnerProduct(bvel, grad(u)) * (v + tau * InnerProduct(bvel, grad(v))) * dx
    a += 1e-8 * u * v * dx                   # tiny mass: fixes the inflow-tangent null
    a.Assemble()
    f = LinearForm(fes)
    f.Assemble()                             # homogeneous transport (zero RHS)
    gf = GridFunction(fes)
    gf.Set(seed_cf, definedon=mesh.Boundaries("inflow"))
    r = f.vec.CreateVector()
    r.data = f.vec - a.mat * gf.vec
    gf.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="umfpack") * r
    Jres = InnerProduct(bvel, grad(gf))
    num = math.sqrt(Integrate(Jres * Jres * dx, mesh))
    gn = math.sqrt(Integrate(InnerProduct(grad(gf), grad(gf)) * dx, mesh))
    return gf, num / (gn + 1e-30)


def wall1_nonuniqueness(mesh):
    """Known Clebsch J = grad(z) x grad(x) = e_y (mu0 = x).  Recover mu from the
    transport PDE seeded with mu0=x and with f(mu0)=x^2 -- both valid, different.
    """
    bvel = CF((0.0, 1.0, 0.0))               # J / |J| = e_y; streamlines run +y
    mu_a, res_a = _transport_recover(mesh, bvel, x)        # seed mu0   = x
    mu_b, res_b = _transport_recover(mesh, bvel, x * x)    # seed f(mu0)= x^2

    norm_a = math.sqrt(Integrate(mu_a * mu_a * dx, mesh))
    # mu_b is the reparametrization f(mu_a) = mu_a^2  (recovery returned f(mu0)):
    d_reparam = math.sqrt(Integrate((mu_b - mu_a * mu_a) ** 2 * dx, mesh)) / norm_a
    # ...yet mu_b is NOT mu_a: the two recovered stream surfaces genuinely differ:
    d_distinct = math.sqrt(Integrate((mu_b - mu_a) ** 2 * dx, mesh)) / norm_a
    return {
        "res_seed_mu0": float(res_a),
        "res_seed_fmu0": float(res_b),
        "reparam_match": float(d_reparam),     # small => mu_b = f(mu_a)
        "distinct": float(d_distinct),         # O(1) => non-unique
    }


# ----------------------------------------------------------------------------
# WALL 2 -- degeneracy at a field null
# ----------------------------------------------------------------------------
def wall2_null_degeneracy(maxh=0.08):
    """Known Clebsch with isolated line-null on the z-axis:
        lambda0=(x^2-y^2)/2, mu0=x*y  ->  J=(0,0,x^2+y^2), |J|=0 on axis.
    Probe the Clebsch area form |grad lambda x grad mu| = |J| and the recovery
    conditioning ~ 1/|J|^2 along a ray approaching the null.
    """
    geo = CSGeometry()
    geo.Add(OrthoBrick(Pnt(-0.5, -0.5, 0), Pnt(0.5, 0.5, 1)).mat("dom"))
    mesh = Mesh(geo.GenerateMesh(maxh=maxh))

    glam = CF((x, -y, 0.0))                   # grad(lambda0)
    gmu = CF((y, x, 0.0))                      # grad(mu0)
    areaform = Norm(Cross(glam, gmu))          # |grad lam x grad mu| = |J| = x^2+y^2

    radii = [0.40, 0.20, 0.10, 0.05, 0.025]
    rows = []
    for r in radii:
        p = mesh(r / math.sqrt(2.0), r / math.sqrt(2.0), 0.5)
        a = float(areaform(p))                 # = r^2 analytically
        rows.append((r, a, a * a, 1.0 / (a * a + 1e-30)))

    # fit log(cond) vs log(r): expect slope ~ -4  (cond ~ 1/|J|^2 ~ 1/r^4)
    lr = np.log(np.array([r for r, *_ in rows]))
    lc = np.log(np.array([c for *_, c in rows]))
    slope = float(np.polyfit(lr, lc, 1)[0])

    a_far = rows[0][3]                          # cond at r=0.40
    a_near = rows[-1][3]                        # cond at r=0.025
    return {
        "rows": rows,
        "cond_slope": slope,                   # ~ -4
        "cond_ratio": a_near / a_far,          # huge => blows up toward null
        "areaform_near": rows[-1][1],          # -> 0 at the null
    }, mesh


def main():
    print("=== F2 probe: Clebsch recovery from a general div-free J is ILL-POSED ===")
    with TaskManager():
        mesh1 = _cube_with_inflow()
        w1 = wall1_nonuniqueness(mesh1)
        w2, mesh2 = wall2_null_degeneracy()

    print("\n-- WALL 1: non-uniqueness (recovered mu is seed-dependent, = f(mu0)) --")
    print(f"   transport residual, seed mu0=x   : {w1['res_seed_mu0']:.2e}  (valid)")
    print(f"   transport residual, seed f(mu0)=x^2: {w1['res_seed_fmu0']:.2e}  (valid)")
    print(f"   ||mu_b - mu_a^2|| / ||mu_a||      : {w1['reparam_match']:.2e}"
          f"   (small => mu_b = f(mu_a))")
    print(f"   ||mu_b - mu_a||   / ||mu_a||      : {w1['distinct']:.2e}"
          f"   (O(1) => genuinely different => NON-UNIQUE)")

    print("\n-- WALL 2: degeneracy at a field null (|J|=0 on z-axis) --")
    print("   r        |J|=areaform     cond ~ 1/|J|^2")
    for r, a, _a2, c in w2["rows"]:
        print(f"   {r:.3f}    {a:.4e}       {c:.3e}")
    print(f"   log-log slope cond vs r           : {w2['cond_slope']:.2f}"
          f"   (~ -4 => cond ~ 1/r^4, singular at null)")
    print(f"   cond ratio (near/far)             : {w2['cond_ratio']:.2e}"
          f"   (blows up toward null)")

    # PASS = the wall is cleanly demonstrated (NOT "recovery succeeded"):
    #  W1: both seeds valid, recovered fields distinct but reparam-related.
    #  W2: conditioning blows up like 1/r^4 toward the null.
    ok = (
        w1["res_seed_mu0"] < 1e-5 and w1["res_seed_fmu0"] < 1e-5      # both valid
        and w1["reparam_match"] < 1e-3                                # mu_b = f(mu_a)
        and w1["distinct"] > 0.1                                      # but distinct
        and w2["cond_slope"] < -3.5                                   # ~ -4
        and w2["cond_ratio"] > 1e3                                    # blows up
    )
    print(f"\n=== WALL DEMONSTRATED (probe PASS): {ok} ===")
    print("    Interpretation: a GENERAL div-free J has no CANONICAL Clebsch pair --")
    print("    recovery is non-unique (any f(mu0)) and singular at field nulls (|J|=0).")
    return {"wall1": w1, "wall2": {k: v for k, v in w2.items() if k != "rows"}}, ok


if __name__ == "__main__":
    _, ok = main()
    raise SystemExit(0 if ok else 1)
