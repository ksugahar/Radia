r"""Air-gap element (AGE) -- NGSolve coupling (step 2).

Wires the analytic gap DtN (airgap_element.annular_dtn_matrix) into an NGSolve FE
solve as a boundary (Robin) operator: the air gap ri<r<ro is NEVER meshed -- only the
stator-side annulus ro<r<Rext is FE, and the gap couples through its harmonic
Dirichlet-to-Neumann action on the gap-facing boundary.

Gate: a single harmonic n, rotor-ring phasor A_in at ri, Dirichlet A_ext*cos(n theta)
at Rext.  The whole ri<r<Rext region is current-free Laplace, so the exact field is
A(r) = alpha r^n + beta r^-n (full-annulus coeffs).  The AGE-coupled FE solution on
ro<r<Rext must reproduce it -- i.e. the un-meshed gap is represented EXACTLY by the
DtN border.  Weak form:

    a(A,v) = (grad A, grad v)_Omega + M22 <A, v>_inner
    L(v)   = -M21 A_in <cos(n theta), v>_inner          (Omega = ro<r<Rext)

with (M21, M22) the stator-ring row of the gap DtN matrix.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from ngsolve import (Mesh, H1, GridFunction, BilinearForm, LinearForm, grad, dx, ds,
                     x, y, atan2, cos, CoefficientFunction, TaskManager)
from netgen.geom2d import SplineGeometry

from radia_mcp.radia_ngsolve.airgap_element import (annular_dtn_matrix,
                                                    annular_harmonic_coeffs)

RI, RO, REXT = 0.10, 0.11, 0.20     # rotor ring / gap-stator ring / FE outer [m]
N = 2
A_IN, A_EXT = 1.0, 0.3              # rotor-ring & outer harmonic coefficients


def annulus_mesh(maxh=0.006):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer")
    geo.AddCircle((0, 0), RO, leftdomain=0, rightdomain=1, bc="inner")
    return Mesh(geo.GenerateMesh(maxh=maxh))


def test_age_dtn_border_reproduces_full_annulus():
    mesh = annulus_mesh()
    mesh.Curve(5)
    cosn = cos(N * atan2(y, x))
    (_, _), (M21, M22) = annular_dtn_matrix(N, RI, RO)        # gap ri->ro stator-ring row

    fes = H1(mesh, order=5, dirichlet="outer")
    A, v = fes.TnT()
    ds_in = ds(definedon=mesh.Boundaries("inner"))
    a = BilinearForm(fes)
    a += grad(A) * grad(v) * dx
    a += M22 * A * v * ds_in                                  # AGE Robin self-term
    f = LinearForm(fes)
    f += (-M21 * A_IN) * cosn * v * ds_in                     # rotor-ring coupling source
    with TaskManager():
        a.Assemble(); f.Assemble()
        gfu = GridFunction(fes)
        gfu.Set(A_EXT * cosn, definedon=mesh.Boundaries("outer"))
        r = f.vec - a.mat * gfu.vec
        gfu.vec.data += a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky") * r

    # exact: full annulus ri<r<Rext, A(ri)=A_in, A(Rext)=A_ext (gap represented exactly)
    alpha, beta = annular_harmonic_coeffs(N, RI, REXT, A_IN, A_EXT)

    diffs, refs = [], []
    for frac in (0.05, 0.4, 0.95):
        rr = RO + frac * (REXT - RO)
        for th in (0.2, 1.0, 2.3, -1.4):
            fe = gfu(mesh(rr * math.cos(th), rr * math.sin(th)))
            exact = (alpha * rr ** N + beta * rr ** -N) * math.cos(N * th)
            diffs.append(abs(fe - exact)); refs.append(abs(exact))
    rel = max(diffs) / max(refs)
    # the un-meshed gap should also hand back the right ring value at ro
    a_out_fe = gfu(mesh(RO * math.cos(0.0), RO * math.sin(0.0)))
    a_out_exact = alpha * RO ** N + beta * RO ** -N
    print(f"AGE/NGSolve: max rel err vs full-annulus analytic = {rel:.2e}  "
          f"(A(ro): FE {a_out_fe:.5f} vs exact {a_out_exact:.5f})")
    assert rel < 5e-3, f"AGE-coupled FE field off analytic by {rel:.2e}"
    assert math.isclose(a_out_fe, a_out_exact, rel_tol=5e-3)


def main():
    test_age_dtn_border_reproduces_full_annulus()
    print("[OK] AGE step 2: the analytic gap DtN, wired as an NGSolve boundary operator, "
          "reproduces the full-annulus field with the gap UN-MESHED. Mesh-free gap coupling works.")


if __name__ == "__main__":
    main()
