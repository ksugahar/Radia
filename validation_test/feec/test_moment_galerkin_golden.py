"""Golden lock for radia.moment_galerkin -- the SYMMETRIC moment-Galerkin MMMM demag solver.

The moment-Galerkin demag operator N = B^T G B (moment basis B = M.n face charges, G = the exact analytic
charge-Gram H-matrix, the same C++ kernel HDiv-VIM ships) is SYMMETRIC by Green-kernel reciprocity, unlike
collocation MMMM (||A - A^T|| ~ 0.07-1.6).  Symmetric -> loop modes field-null by construction -> loop-free /
mu_r-independent convergence with no loop-star.

These tests lock the de-risk-validated dipole-level solver (3 DOF/hex):
  (a) cube operator demag factor d ~ 1/3;
  (b) single-cube self-consistency M_z == chi/(1+d chi) H0 across mu_r (machine precision);
  (c) N = B^T G B is SYMMETRIC to machine precision (the core loop-free claim) on a multi-hex block;
  (d) a 2-hex bar's iron STRAY field matches rad.Solve (full 6-DOF yano-MSC) at a far probe (dipole vs 6-DOF).

Self-contained (mesh-less hex vertex lists + the existing C++ charge-Gram + mass-Riesz CG), fast, no NGSolve.
"""
import numpy as np
import pytest
import radia as rad
import radia.moment_galerkin as mg

MU0 = 4e-7 * np.pi
H0 = 1000.0
D = 0.02


def _hexv(x0, x1, y0, y1, z0, z1):
    return np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                     [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]], float)


@pytest.fixture(autouse=True)
def _clean():
    rad.UtiDelAll(); rad.set_demag_backend("auto")
    yield
    rad.set_demag_backend("auto"); rad.UtiDelAll()


def test_cube_demag_factor():
    """(a) the cube operator demag factor is ~1/3 (the residual ~0.4% is the coincident-face Gram quad)."""
    sys = mg.assemble_moment_system([_hexv(0, D, 0, D, 0, D)])
    d = mg.demag_factor(sys, 2)
    assert 0.330 < d < 0.340, f"cube demag factor {d} out of [0.330, 0.340]"


def test_material_input_validation():
    """Invalid soft-iron susceptibility is rejected before the C++ solve sees 1/chi."""
    hexes = [_hexv(0, D, 0, D, 0, D)]
    with pytest.raises(ValueError, match="mu_r > 1|chi > 0"):
        mg.moment_galerkin_demag_solve(hexes, mu_r=1.0, H_ext=(0.0, 0.0, H0))
    sys = mg.assemble_moment_system(hexes)
    with pytest.raises(ValueError, match="chi must be positive"):
        mg.solve_assembled(sys, (0.0, 0.0, H0), 0.0)


@pytest.mark.parametrize("mu_r", [10.0, 100.0, 1000.0])
def test_cube_self_consistency(mu_r):
    """(b) the cube solve reproduces the closed-form M_z = chi/(1 + d chi) H0 to machine precision."""
    chi = mu_r - 1.0
    sys = mg.assemble_moment_system([_hexv(0, D, 0, D, 0, D)])
    d = mg.demag_factor(sys, 2)
    M, _ = mg.solve_assembled(sys, (0.0, 0.0, H0), chi)
    pred = chi / (1.0 + d * chi) * H0
    assert abs(M[0, 2] - pred) <= 1e-6 * abs(pred), f"M_z {M[0,2]} != {pred}"


def test_operator_symmetry():
    """(c) N = B^T G B is symmetric to machine precision (the loop-free claim) on a 2x2x1 block."""
    hexes = [_hexv(i * D, (i + 1) * D, j * D, (j + 1) * D, 0, D) for i in range(2) for j in range(2)]
    sys = mg.assemble_moment_system(hexes)
    G, B = sys["G"], sys["B"]
    ndof = B.shape[1]
    N = np.empty((ndof, ndof))
    Bc = B.tocsr()
    for j in range(ndof):
        ej = np.zeros(ndof); ej[j] = 1.0
        c = np.asarray(Bc @ ej)
        N[:, j] = Bc.T @ np.asarray(G.matvec(c.tolist()), float)
    asym = np.linalg.norm(N - N.T) / np.linalg.norm(N)
    assert asym < 1e-9, f"moment-Galerkin N not symmetric: ||N-N^T||/||N|| = {asym:.2e}"


def test_bar_stray_field_vs_radsolve():
    """(d) a 2-hex bar's iron stray field matches rad.Solve (full 6-DOF yano-MSC) at a far probe to ~1%."""
    bar = [_hexv(0, D, 0, D, 0, D), _hexv(D, 2 * D, 0, D, 0, D)]
    mu_r = 1000.0
    probe = [0.0, 0.0, 5 * D]
    out = mg.moment_galerkin_demag_solve(bar, mu_r=mu_r, H_ext=(0.0, 0.0, H0))
    M = out["M"]
    rad.UtiDelAll()
    objs = [rad.ObjHexahedron([list(v) for v in V], list(M[e])) for e, V in enumerate(bar)]
    B_mg = np.array(rad.Fld(rad.ObjCnt(objs), "b", probe))
    rad.UtiDelAll(); rad.set_demag_backend("yano")
    objs2 = [rad.ObjHexahedron([list(v) for v in V], [0, 0, 0]) for V in bar]
    for h in objs2:
        rad.MatApl(h, rad.MatLin(mu_r))
    cont = rad.ObjCnt(objs2 + [rad.ObjBckg(lambda p: [0, 0, MU0 * H0])])
    rad.Solve(cont, 1e-6, 2000, 0)
    B_rad = np.array(rad.Fld(cont, "b", probe)) - np.array([0.0, 0.0, MU0 * H0])
    rel = np.linalg.norm(B_mg - B_rad) / max(np.linalg.norm(B_rad), 1e-30)
    assert rel < 1e-2, f"bar stray field rel diff {rel:.2e} vs rad.Solve too large"
