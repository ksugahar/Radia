r"""Weak-form pullback kata: the 2-D Laplace weak form in differential forms,
pulled into the (A, phi) potential chart -> the hodograph weight W = diag(mu, 1/mu).

This is the "kata" companion to bidirectional_coordinate_transform_2d.py.  Where
that file proves a 2-D CONFORMAL map is weight-free (W = I) for a GENERIC map and
verifies it with a deformed-mesh NGSolve assembly, this file works the SAME
exterior-calculus pullback for the SPECIFIC magnetostatic conjugate map -- the
hodograph -- and shows the material mu survives as the ANISOTROPIC weight
W = diag(mu, 1/mu), i.e. exactly the Tampere bidirectional governing equation (8).
It is sympy-only (fast, no FE solve): a teaching + golden artifact for the
"weak form in differential forms" practice.

THE KATA (the 6 steps; scalar Laplace in forms)
-----------------------------------------------
Unknown 0-form u, material mu (the Hodge star carries the material):
  1. strong:   d(mu * star d u) = 0
  2. test:     0-form v
  3. Leibniz:  d(v mu star du) = dv ^ mu star du + v d(mu star du)
  4. Stokes:   INT_Omega d(v mu star du) = INT_dOmega v mu star du
  5. read off: a(u,v) = INT_Omega mu dv ^ star du = INT_Omega mu grad u . grad v
               natural BC = INT_dOmega v mu d_n u   (Neumann flux)
  6. element:  0-form -> H1 (nodal)
The single fact that makes this ONE line is  alpha ^ star beta = <alpha,beta> vol
(the wedge-with-Hodge IS the L2 inner product = the FE mass/stiffness; d is the
metric-free part, the Hodge star carries the material).

LV4 -- pull the weak form into the (A, phi) potential chart
----------------------------------------------------------
The weak form INT dv ^ mu star du is coordinate-free, so a chart change moves
ONLY the Hodge star (the metric), as the weight W = |det J| (J^T J)^{-1}.  For the
magnetostatic conjugate pair (A = flux function, phi = scalar potential) the
equipotentials phi=const and flux lines A=const are orthogonal (Hodge star), so
(A, phi) is an orthogonal chart.  With field magnitude q = |grad phi| = |H| and
|grad A| = |B| = mu q, the physical line element is

    ds^2 = dA^2 / (mu q)^2 + dphi^2 / q^2   =>   g = diag(1/(mu q)^2, 1/q^2),

so the FE weight  W = sqrt(det g) g^{-1} = diag(mu, 1/mu) -- INDEPENDENT of the
field magnitude q (the conformal factor cancels).  That q-independence is WHY the
hodograph linearises: for nonlinear mu(q) the weight is W = diag(mu(q), 1/mu(q)),
a KNOWN function of the coordinate (Chaplygin).

So the (A, phi) Laplace is, by the SAME 6 steps,
    d_A(mu d_A x) + d_phi((1/mu) d_phi x) = 0   <=>   div_(A,phi)(W grad x) = 0,
which is the Tampere bidirectional equation (8); the unknown is now the GEOMETRY
x(A,phi), y(A,phi) -- independent and dependent variables have SWAPPED (the
hodograph), valid where the chart Jacobian d(x,y)/d(A,phi) != 0.

VERIFIED HERE (symbolic, sympy)
-------------------------------
  * line-element derivation: g = diag(1/(mu q)^2, 1/q^2) -> W = diag(mu, 1/mu),
    with the field magnitude q cancelling exactly.
  * conjugate-map pullback for the conformal maps F(z) = z, z^2, z^3, 1/z
    (forward chart map (A, phi) = (Re F, Im F / mu)):  W = diag(mu, 1/mu) EXACTLY
    in every case -- the conformal factor cancels regardless of the map.

run:  python weakform_pullback_kata.py
"""

from _validation_output import validation_output
import os

import sympy as sp


def symbolic_W_from_line_element():
    """W from the orthogonal (A, phi) line element ds^2 = dA^2/(mu q)^2 + dphi^2/q^2.

    Returns the diagonal of W = sqrt(det g) g^{-1} and whether it equals
    diag(mu, 1/mu) with the field magnitude q cancelled."""
    mu, q = sp.symbols("mu q", positive=True)
    g = sp.diag(1 / (mu * q) ** 2, 1 / q ** 2)
    W = sp.simplify(sp.sqrt(g.det()) * g.inv())
    target = sp.diag(mu, 1 / mu)
    return {
        "W_diag": [str(W[0, 0]), str(W[1, 1])],
        "q_cancels": bool(sp.simplify(W - target) == sp.zeros(2, 2)),
    }


def _W_for_conjugate_map(F):
    """W = |det J| (J^T J)^{-1} for the chart map (A, phi) = (Re F(z), Im F(z)/mu).

    z = x + i y;  J = d(x,y)/d(A,phi) = J_fwd^{-1}, so
    (J^T J)^{-1} = J_fwd J_fwd^T and |det J| = 1/|det J_fwd|, where
    J_fwd = d(A,phi)/d(x,y)."""
    x, y = sp.symbols("x y", real=True)
    mu = sp.symbols("mu", positive=True)
    z = x + sp.I * y
    Fz = F(z)
    A = sp.re(Fz)
    phi = sp.im(Fz) / mu
    Jf = sp.Matrix([[sp.diff(A, x), sp.diff(A, y)],
                    [sp.diff(phi, x), sp.diff(phi, y)]])
    detJf = sp.simplify(Jf.det())
    W = sp.simplify((1 / sp.Abs(detJf)) * (Jf * Jf.T))
    return W, mu


MAPS = [("z", lambda z: z),
        ("z^2", lambda z: z ** 2),
        ("z^3", lambda z: z ** 3),
        ("1/z", lambda z: 1 / z)]


def run():
    """Symbolic kata verification; returns a dict for the golden test."""
    le = symbolic_W_from_line_element()
    cases = []
    for name, F in MAPS:
        W, mu = _W_for_conjugate_map(F)
        target = sp.Matrix([[mu, 0], [0, 1 / mu]])
        cases.append({
            "F": name,
            "W_diag": [str(W[0, 0]), str(W[1, 1])],
            "W_offdiag": [str(W[0, 1]), str(W[1, 0])],
            "is_diag_mu_invmu": bool(sp.simplify(W - target) == sp.zeros(2, 2)),
        })
    return {
        "line_element": le,
        "cases": cases,
        "all_diag_mu_invmu": all(c["is_diag_mu_invmu"] for c in cases),
    }


def main():
    import json
    print("=" * 72)
    print("Weak-form pullback kata: 2-D Laplace in forms -> hodograph W = diag(mu, 1/mu)")
    print("=" * 72)
    out = run()
    le = out["line_element"]
    print("\nLv4 line-element derivation  g = diag(1/(mu q)^2, 1/q^2):")
    print("  W = diag(%s, %s)   q cancels: %s"
          % (le["W_diag"][0], le["W_diag"][1], le["q_cancels"]))
    print("\nconjugate-map pullback  (A, phi) = (Re F, Im F / mu):")
    print("  %-6s %-26s %-12s %s" % ("F(z)", "W diag", "W offdiag", "== diag(mu,1/mu)?"))
    for c in out["cases"]:
        print("  %-6s diag(%s, %s)%s(%s,%s)   %s"
              % (c["F"], c["W_diag"][0], c["W_diag"][1],
                 " " * max(1, 18 - len(c["W_diag"][0]) - len(c["W_diag"][1])),
                 c["W_offdiag"][0], c["W_offdiag"][1], c["is_diag_mu_invmu"]))
    print("\nall conformal maps give W = diag(mu, 1/mu):", out["all_diag_mu_invmu"])

    here = os.path.dirname(__file__)
    with validation_output("weakform_pullback_kata.json").open("w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2)
    print("\nsaved weakform_pullback_kata.json")


if __name__ == "__main__":
    main()
