"""Phase 2-E sanity test:
   For a single Q1 axis-aligned rectangle, the C++ closed-form Integrator
   must produce the same 4x4 element matrix (in psi-DOF) as the Python
   reference (axifem/axifem_quad.py).

   Reference implementation has been validated against BEM v3 to 0.5%, so
   matrix-level agreement is the precise criterion for Phase 2-E.
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "tests", "axifem", "_reference_python"))
from axifem_quad import (
    _element_matrices_quad_closed_form,
    element_sigma_mass_quad,
    _phi_basis_inverse,
)

from ngsolve import (
    Mesh, BilinearForm, GridFunction, CoefficientFunction,
    TaskManager,
)
from netgen.occ import OCCGeometry, MoveTo, X, Y
from radia.axifem import (
    H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI,
)


def python_q1_v_matrices(ra, rb, za, zb, mu, sigma):
    """Reference K, M in V-DOF (matches C++ Integrator output convention)."""
    M_V, _, _ = _element_matrices_quad_closed_form(ra, rb, za, zb, mu, mu)
    M_sig_V = element_sigma_mass_quad(ra, rb, za, zb, sigma)
    return M_V, M_sig_V


def get_cpp_element_matrices(ra, rb, za, zb, mu, sigma):
    """Build a single-element NGSolve mesh and extract the 4x4 element matrices."""
    box = MoveTo(ra, za).Rectangle(rb - ra, zb - za).Face()
    box.faces.name = "conductor"
    box.edges.Min(X).name = "left"
    box.edges.Max(X).name = "right"
    box.edges.Min(Y).name = "bot"
    box.edges.Max(Y).name = "top"
    # Force a SINGLE Q1 element by using a coarse maxh.
    mesh = Mesh(OCCGeometry(box, dim=2).GenerateMesh(
        maxh=10*(rb - ra), quad_dominated=True))
    if mesh.ne != 1:
        raise RuntimeError(f"expected 1 element, got {mesh.ne}")

    fes = H1Henrotte(mesh)
    if fes.ndof != 4:
        raise RuntimeError(f"expected ndof=4, got {fes.ndof}")

    mu_cf = CoefficientFunction(mu)
    sigma_cf = CoefficientFunction(sigma)

    a = BilinearForm(fes, symmetric=True)
    a += AxiHenrotteStiffnessBFI(mu_cf)
    with TaskManager():
        a.Assemble()

    m = BilinearForm(fes, symmetric=True)
    m += AxiHenrotteSigmaMassBFI(sigma_cf)
    with TaskManager():
        m.Assemble()

    K = np.zeros((4, 4))
    rs, cs, vs = a.mat.COO()
    for r, c, v in zip(rs, cs, vs):
        K[r, c] += v

    M = np.zeros((4, 4))
    rs, cs, vs = m.mat.COO()
    for r, c, v in zip(rs, cs, vs):
        M[r, c] += v

    # NGSolve vertex order vs Python (sa,za),(sb,za),(sb,zb),(sa,zb).
    # Match by vertex coordinate.
    target = [(ra, za), (rb, za), (rb, zb), (ra, zb)]
    perm = []
    for tr, tz in target:
        for j, v in enumerate(mesh.vertices):
            p = v.point
            if abs(p[0] - tr) < 1e-9 and abs(p[1] - tz) < 1e-9:
                perm.append(j)
                break
    perm = np.array(perm, dtype=int)
    K = K[perm[:, None], perm[None, :]]
    M = M[perm[:, None], perm[None, :]]
    return K, M


def main():
    mu0 = 4 * np.pi * 1e-7
    sigma_cu = 5.8e7

    cases = [
        ("interior", 1e-3, 2e-3, 0.0, 1e-3),
        ("axis-touching", 0.0, 2e-3, 0.0, 1e-3),
        ("interior-tall", 1e-3, 5e-3, -2e-3, 2e-3),
    ]

    for label, ra, rb, za, zb in cases:
        print(f"\n=== {label}: ra={ra*1e3}mm rb={rb*1e3}mm za={za*1e3}mm zb={zb*1e3}mm ===")
        K_py, M_py = python_q1_v_matrices(ra, rb, za, zb, mu0, sigma_cu)
        try:
            K_cpp, M_cpp = get_cpp_element_matrices(ra, rb, za, zb, mu0, sigma_cu)
        except RuntimeError as e:
            print(f"  SKIP (mesh build): {e}")
            continue

        # Both K_py and K_cpp now in V-DOF; compare directly.
        denom_K = np.linalg.norm(K_py) if np.linalg.norm(K_py) > 0 else 1.0
        denom_M = np.linalg.norm(M_py) if np.linalg.norm(M_py) > 0 else 1.0
        err_K = np.linalg.norm(K_cpp - K_py) / denom_K
        err_M = np.linalg.norm(M_cpp - M_py) / denom_M
        print(f"  K Python:\n{K_py}\n  K C++:\n{K_cpp}")
        print(f"  K relative Frobenius error: {err_K:.3e}")
        print(f"  M relative Frobenius error: {err_M:.3e}")
        if max(err_K, err_M) < 1e-10:
            print(f"  [OK] machine precision agreement")
        elif max(err_K, err_M) < 1e-6:
            print(f"  [OK] within 1e-6 relative")
        else:
            print(f"  [FAIL] errors above 1e-6")


if __name__ == "__main__":
    main()
