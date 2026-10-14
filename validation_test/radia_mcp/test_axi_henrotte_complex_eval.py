"""COMPLEX pointwise field evaluation on the AxiHenrotte (H1Henrotte) space.

Before the C++ fix, evaluating a COMPLEX H1Henrotte GridFunction (value OR
gradient) fell back to the NGSolve base-class DifferentialOperator stub and
returned wrong values -- only the REAL CalcMatrix overloads of AxiHenrotteDiffOpId
/ AxiHenrotteDiffOpGradient were implemented.  This blocked any field-based
post-processing of solve_axi_eddy_harmonic (|A|, |B|), and hence the nonlinear
mu(|B|) Picard layer.  src/ext/axifem/axi_henrotte_diffop.hpp now also provides
the BareSliceMatrix<Complex,ColMajor> CalcMatrix overloads (shape/dshape are
real; only the DOF coefficients are complex), so complex eval mirrors the real
path exactly.

Test idea (isolates the EVAL, independent of .Set / projection): build the SAME
nodal field on a real and on a complex H1Henrotte space, with the complex DOFs
set to ``c`` times the real DOFs.  Then at every off-axis point the complex
evaluation MUST equal ``c`` times the real evaluation -- for the value (Id) and
for the gradient (Gradient).  A broken (base-stub) complex eval fails this.
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
from netgen.occ import MoveTo, OCCGeometry, X, Y
from ngsolve import GridFunction, Mesh, TaskManager, grad
from ngsolve import x as r_cf
from ngsolve import y as z_cf
from radia.axifem import H1Henrotte

C = 2.0 + 3.0j
PTS = [(0.30, 0.40), (0.50, 0.55), (0.72, 0.20), (0.88, 0.83)]


def build_mesh():
    rect = MoveTo(0, 0).Rectangle(1.0, 1.0).Face()
    rect.edges.Min(X).name = "axis"
    rect.edges.Max(X).name = "outer"
    rect.edges.Min(Y).name = "outer"
    rect.edges.Max(Y).name = "outer"
    return Mesh(OCCGeometry(rect, dim=2).GenerateMesh(maxh=0.12))


def main():
    mesh = build_mesh()

    # A_phi = r^2 (1 + z) is inside the Q1 Henrotte basis
    # {1, r^2, z, r^2 z} and has two nonzero gradient components.  Using
    # r^2 alone makes a componentwise relative check divide roundoff in the
    # analytically zero z derivative by an almost-zero reference value.
    f = r_cf * r_cf * (1.0 + z_cf)

    fes_r = H1Henrotte(mesh, order=1, dirichlet="")
    gfr = GridFunction(fes_r)
    gfr.Set(f)

    fes_c = H1Henrotte(mesh, order=1, complex=True, dirichlet="")
    gfc = GridFunction(fes_c)
    # complex DOFs := C * real DOFs  -> the complex field is exactly C * (real field)
    gfc.vec.FV().NumPy()[:] = C * np.asarray(gfr.vec.FV().NumPy())

    print(f"  C = {C}")
    print("  --- VALUE (Id) eval: gfc(pt) vs C*gfr(pt) ---")
    max_v = 0.0
    for (rr, zz) in PTS:
        vr = gfr(mesh(rr, zz))
        vc = gfc(mesh(rr, zz))
        err = abs(vc - C * vr) / (abs(C * vr) + 1e-30)
        max_v = max(max_v, err)
        print(f"    ({rr:.2f},{zz:.2f}): gfr={vr:+.5e}  gfc={vc.real:+.4e}{vc.imag:+.4e}j"
              f"  rel={err:.2e}")
    assert max_v < 1e-9, f"COMPLEX value eval != C * real eval (max rel {max_v:.2e})"

    print("  --- GRADIENT eval: grad(gfc)(pt) vs C*grad(gfr)(pt) ---")
    max_g = 0.0
    for (rr, zz) in PTS:
        gr = grad(gfr)(mesh(rr, zz))     # (dA/dr, dA/dz) real
        gc = grad(gfc)(mesh(rr, zz))     # complex
        expected = np.asarray(C * np.asarray(gr), dtype=np.complex128)
        actual = np.asarray(gc, dtype=np.complex128)
        err = np.linalg.norm(actual - expected) / max(np.linalg.norm(expected), 1e-30)
        max_g = max(max_g, float(err))
        print(f"    ({rr:.2f},{zz:.2f}): grad_r=({gr[0]:+.4e},{gr[1]:+.4e})"
              f"  |grad_c-C*grad_r|rel<= {max_g:.2e}")
    assert max_g < 1e-9, f"COMPLEX gradient eval != C * real gradient (max rel {max_g:.2e})"

    # ensure the comparison is not vacuous: gfr must be a non-trivial field (not
    # all-zero), otherwise gfc == C*gfr would hold trivially as 0 == 0.
    assert abs(gfr(mesh(0.5, 0.55))) > 1e-2, "real field is trivially ~zero"

    print(f"\n[OK] COMPLEX H1Henrotte field eval works: value (Id) and gradient "
          f"(Gradient)\n     both equal C * the real eval to machine precision "
          f"(max rel value {max_v:.1e}, grad {max_g:.1e}).\n     The base-class "
          f"stub is no longer hit -- complex post-processing of "
          f"solve_axi_eddy_harmonic\n     (|A|, |B|) is now reliable, unblocking "
          f"the nonlinear mu(|B|) Picard layer.")


if __name__ == "__main__":
    with TaskManager():
        main()
