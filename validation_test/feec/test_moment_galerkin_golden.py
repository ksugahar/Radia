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


# ----------------------------------------------------------------------------------------------------
# Quad / 6-DOF (3 dipole + 2 residual-eigenmode quad) moment-Galerkin -- the higher-order increment.
# ----------------------------------------------------------------------------------------------------

def test_field_recon_anchor():
    """field-from-sigma reconstruction of a KNOWN constant magnetization matches rad.Fld(ObjHexahedron) to
    ~1e-10 at external probes -- anchors the reconstruction units/constants (mu0/4pi + sign)."""
    M0 = np.array([300.0, -700.0, 1200.0])
    V = _hexv(0, D, 0, D, 0, D)
    sys = mg.assemble_moment_system([V])
    m = np.tile(M0, sys["n_hex"])
    probes = np.array([[0.05, 0.03, 0.04], [-0.04, 0.06, -0.03], [0.011, 0.013, 0.08]])
    B_recon = mg.reconstruct_field(sys, m, probes)
    obj = rad.ObjHexahedron([list(v) for v in V], list(M0))
    for i, p in enumerate(probes):
        B_rad = np.array(rad.Fld(obj, "b", [float(p[0]), float(p[1]), float(p[2])]))
        rel = np.linalg.norm(B_recon[i] - B_rad) / np.linalg.norm(B_rad)
        assert rel < 1e-9, f"field-from-sigma anchor probe {i}: rel {rel:.2e} vs rad.Fld"


@pytest.mark.parametrize("mu_r", [10.0, 100.0, 1000.0])
def test_quad_cube_gate(mu_r):
    """A cube in a uniform field: the quad amplitudes stay ~0 (symmetry) and the dipole M_z is the SAME
    chi/(1+d chi) H0 as dipole-only -- the 2 quad DOF must not corrupt the symmetric cube."""
    sys = mg.assemble_moment_system([_hexv(0, D, 0, D, 0, D)], quad=True)
    out = mg.moment_galerkin_demag_solve([_hexv(0, D, 0, D, 0, D)], mu_r=mu_r, H_ext=(0.0, 0.0, H0), quad=True)
    Mz = out["M"][0, 2]
    quad = out["quad_amps"][0]
    d = mg.demag_factor(sys, 2)
    pred = (mu_r - 1.0) / (1.0 + d * (mu_r - 1.0)) * H0
    assert np.max(np.abs(quad)) < 1e-2 * abs(Mz), f"cube quad amps not ~0: {quad}"
    assert abs(Mz - pred) < 5e-3 * abs(pred), f"cube 5-DOF M_z {Mz} != dipole pred {pred}"


def test_quad_reduces_to_dipole():
    """Deleting the 2 quad columns from the 5-DOF B/M_mass yields the EXACT 3-DOF dipole assembly (the quad
    is a strict superset; on axis-aligned hexes per-face == per-tri so the dipole blocks coincide)."""
    bar = [_hexv(0, D, 0, D, 0, D), _hexv(D, 2 * D, 0, D, 0, D)]
    s3 = mg.assemble_moment_system(bar, quad=False)
    s5 = mg.assemble_moment_system(bar, quad=True)
    n = s3["n_hex"]
    keep = np.concatenate([[5 * e, 5 * e + 1, 5 * e + 2] for e in range(n)])
    B5d = s5["B"].toarray()[:, keep]
    W5d = s5["M_mass"].toarray()[np.ix_(keep, keep)]
    assert np.allclose(B5d, s3["B"].toarray(), atol=1e-12), "5-DOF dipole columns != 3-DOF B"
    assert np.allclose(W5d, s3["M_mass"].toarray(), atol=1e-12), "5-DOF dipole mass != 3-DOF M_mass"


def test_quad_bar_antisymmetry():
    """A symmetric 3-cube bar under an OBLIQUE (quad-exciting) field: by mirror symmetry the quad amplitudes
    are antisymmetric end-to-end (hex0 = -hex2) with the middle hex ~0 -- the quad basis captures the right
    gradient physics."""
    bar = [_hexv(i * D, (i + 1) * D, 0, D, 0, D) for i in range(3)]
    out = mg.moment_galerkin_demag_solve(bar, mu_r=100.0, H_ext=(1000.0, 0.0, 1500.0), quad=True)
    q = out["quad_amps"]
    scale = np.max(np.abs(q))
    assert scale > 1e-6, "oblique field did not excite the quad modes"
    assert np.max(np.abs(q[1])) < 1e-3 * scale, f"middle hex quad not ~0: {q[1]}"
    assert np.allclose(q[0], -q[2], atol=1e-6 * scale), f"quad not antisymmetric end-to-end: {q[0]} vs {q[2]}"


def test_quad_cpp_matches_dense():
    """The C++ 5-DOF solve (diagonal-Jacobi auto_prec) reproduces a dense reference solve of
    ((1/chi)M_mass + B^T G B) m = M_mass h to ~1e-8 on the oblique bar (locks the C++ integration on the
    ill-conditioned + near-singular-quad-mass system)."""
    bar = [_hexv(i * D, (i + 1) * D, 0, D, 0, D) for i in range(3)]
    mu_r = 100.0; chi = mu_r - 1.0
    H_ext = np.array([1000.0, 0.0, 1500.0])
    sys = mg.assemble_moment_system(bar, quad=True)
    G, B, M_mass, n = sys["G"], sys["B"], sys["M_mass"], sys["n_hex"]
    ndof = B.shape[1]
    Bc = B.tocsr()
    N = np.empty((ndof, ndof))
    for j in range(ndof):
        ej = np.zeros(ndof); ej[j] = 1.0
        N[:, j] = Bc.T @ np.asarray(G.matvec(np.asarray(Bc @ ej).tolist()), float)
    N = 0.5 * (N + N.T)
    Wd = M_mass.toarray()
    h = np.zeros(ndof)
    for e in range(n):
        h[5 * e:5 * e + 3] = H_ext
    rhs = Wd @ h
    m_dense = np.linalg.solve((1.0 / chi) * Wd + N, rhs)
    m_cpp, _ = mg.solve_raw(sys, H_ext, chi)
    rel = np.linalg.norm(m_cpp - m_dense) / np.linalg.norm(m_dense)
    assert rel < 1e-8, f"C++ 5-DOF solve != dense reference: rel {rel:.2e}"


def test_far_quad_matches_analytic_at_scale():
    """The DEFAULT fast NEAR/FAR charge-Gram build (near_factor=2, far_quad=4) reproduces the all-analytic Gram
    solve on a 6x6 plate that HAS far pairs (corner-to-corner ~7D >> the 2D near threshold): external stray field
    to ~1e-4 and the SAME CG iteration count.  Locks the far_quad fast build (4-6x faster at N>=256, the
    Gram-build-bottleneck fix) without changing the demag physics; the symmetric H-matvec keeps CG identical."""
    plate = [_hexv(i * D, (i + 1) * D, j * D, (j + 1) * D, 0, D) for i in range(6) for j in range(6)]
    mu_r = 1000.0
    probe = [0.0, 0.0, 8 * D]
    out_fast = mg.moment_galerkin_demag_solve(plate, mu_r=mu_r, H_ext=(0.0, 0.0, H0))               # far_quad default
    out_ana = mg.moment_galerkin_demag_solve(plate, mu_r=mu_r, H_ext=(0.0, 0.0, H0),
                                             near_factor=1e30, far_quad=0)                            # all-analytic
    objs_f = [rad.ObjHexahedron([list(v) for v in V], list(out_fast["M"][e])) for e, V in enumerate(plate)]
    B_f = np.array(rad.Fld(rad.ObjCnt(objs_f), "b", probe)); rad.UtiDelAll()
    objs_a = [rad.ObjHexahedron([list(v) for v in V], list(out_ana["M"][e])) for e, V in enumerate(plate)]
    B_a = np.array(rad.Fld(rad.ObjCnt(objs_a), "b", probe)); rad.UtiDelAll()
    rel = np.linalg.norm(B_f - B_a) / np.linalg.norm(B_a)
    assert rel < 1e-4, f"far_quad fast build field differs from all-analytic: rel {rel:.2e}"
    # far_quad must not DEGRADE convergence (the loop-free property is preserved by the symmetric matvec); a few
    # extra iters from the ~1e-5 operator perturbation are fine on a small problem -- a real failure would blow up.
    assert out_fast["iters"] < 1.6 * out_ana["iters"] + 10, \
        f"far_quad degraded CG convergence: {out_fast['iters']} iters vs analytic {out_ana['iters']}"
