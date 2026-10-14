"""Golden tests for the projected (Schur-form) mixed Galerkin admittance.

Locks, on the copper cube of validation_test/mixed_galerkin/cube3d:
  * the graded 1-D rule against the closed-form layer integrals,
  * the grid surface block K_ss and b_s against the closed forms,
  * the drive projection against Integrate(phi_n f) (the free-dof-block bug),
  * DC exactness, the wall-band value against the Aitken Foster reference,
    the high-frequency value against the Mellin asymptote,
  * the multi-port form: symmetry, decoupled dipole port, DC Gram.
"""
import cmath
import json
import math
from pathlib import Path

import numpy as np
import pytest

from radia.maglev.mixed_galerkin import schur as SC
from radia.maglev.mixed_galerkin import alpha as A

SIGMA_CU = 5.8e7
MU0 = 4 * math.pi * 1e-7
L = 5e-3
MS = MU0 * SIGMA_CU

# Aitken-extrapolated Foster sums (N = 799) of the exact triple-sine series,
# validation_test/mixed_galerkin/cube3d/04_rank_N_bulk_sweep.py::Y_AITKEN
Y_AITKEN = {1e3: 6.892354, 1e4: 3.431919, 1e5: 1.218962, 1e6: 0.399611}
ROOT = Path(__file__).resolve().parents[1]


def _t(f):
    return cmath.sqrt(2j * math.pi * f * MS)


def test_axis_rule_matches_closed_forms():
    nodes, weights = SC.graded_axis_rule(L)
    assert len(nodes) == 108
    assert nodes.min() > 0.0 and nodes.max() < L
    assert weights.sum() == pytest.approx(L, rel=1e-13)
    for f in (1.0, 1e2, 1e4, 1e6, 1e8, 1e9):
        e0, e2, d2 = SC.check_axis_rule(L, _t(f), nodes, weights)
        assert e0 < 1e-5 and e2 < 1e-5 and d2 < 3e-4, (f, e0, e2, d2)


def test_layer_closed_forms_against_quadrature_at_low_t():
    """box_layer_1d is exact: cross-check with scipy on a smooth case."""
    from scipy.integrate import quad
    t = _t(10.0)
    F0, F2, D2 = SC.box_layer_1d(t, L)

    def part(fun, which):
        re = quad(lambda u: getattr(fun(np.array([u]))[0], which), 0, L)[0]
        return re

    f_re = part(lambda u: SC.layer_values(u, t, L)[0], "real")
    f_im = part(lambda u: SC.layer_values(u, t, L)[0], "imag")
    assert complex(f_re, f_im) == pytest.approx(F0, rel=1e-9)
    f2_re = part(lambda u: SC.layer_values(u, t, L)[0] ** 2, "real")
    f2_im = part(lambda u: SC.layer_values(u, t, L)[0] ** 2, "imag")
    assert complex(f2_re, f2_im) == pytest.approx(F2, rel=1e-9)
    d2_re = part(lambda u: SC.layer_values(u, t, L)[1] ** 2, "real")
    d2_im = part(lambda u: SC.layer_values(u, t, L)[1] ** 2, "imag")
    assert complex(d2_re, d2_im) == pytest.approx(D2, rel=1e-9)


@pytest.fixture(scope="module")
def cube():
    pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")
    from ngsolve import Mesh, TaskManager

    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    for f in box.faces:
        f.name = "outer"
    with TaskManager():
        mesh = Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=L / 5))
        mg = SC.BoxMixedGalerkin(mesh, SIGMA_CU, MU0, n_eigen=12,
                                 rule={"ratio": 2.5})
    return mesh, mg


def test_drive_projection_equals_integrate(cube):
    """The bulk overlaps b_n = int phi_n f dV use the FULL mass matrix."""
    from ngsolve import GridFunction, Integrate, TaskManager, x
    mesh, mg = cube
    with TaskManager():
        lam, vecs, M, free, fes, V = A._dirichlet_eigenmodes(mesh, 6, "outer")
        for cf in (1.0, x - L / 2):
            b = A._project_drive(fes, M, free, vecs, cf, V) * math.sqrt(V)
            gf = GridFunction(fes)
            full = np.zeros(fes.ndof)
            for n in range(vecs.shape[1]):
                full[:] = 0.0
                full[free] = vecs[:, n]
                gf.vec.FV().NumPy()[:] = full
                ref = float(Integrate(gf * cf, mesh, order=6))
                assert b[n] == pytest.approx(ref, abs=1e-12 * abs(b[0]) + 1e-16)
    # and the first overlap is the one of the sine mode, (2/L)^1.5 (2L/pi)^3
    assert abs(mg.B[0, 0]) == pytest.approx((2 / L) ** 1.5 * (2 * L / math.pi) ** 3, rel=2e-3)


def test_surface_block_matches_closed_forms(cube):
    _, mg = cube
    for f in (1e2, 1e4, 1e6):
        s = 2j * math.pi * f
        t = _t(f)
        F0, F2, D2 = SC.box_layer_1d(t, L)
        _, K_ss, b_s, _ = mg.blocks(s)
        K_ref = 3 * D2 * F2 ** 2 + t * t * F2 ** 3
        assert K_ss[0, 0] == pytest.approx(K_ref, rel=2e-3)
        assert b_s[0, 0] == pytest.approx(F0 ** 3, rel=1e-4)


def test_grid_reproduces_mode_normalisation_and_volume(cube):
    _, mg = cube
    assert mg.W3.sum() == pytest.approx(L ** 3, rel=1e-12)
    assert np.sum(mg.W3 * mg.phi[0] ** 2) == pytest.approx(1.0, rel=2e-4)
    assert mg.V == pytest.approx(L ** 3, rel=1e-12)


def test_dc_is_exact_and_band_matches_references(cube):
    _, mg = cube
    Y_DC = SIGMA_CU * L ** 3
    assert mg.Y(0) == pytest.approx(Y_DC, rel=1e-12)
    assert abs(mg.Y(2j * math.pi * 1.0)) == pytest.approx(Y_DC, rel=1e-3)
    # wall band and above: the rank-N mixed Galerkin floor is +0.33 % at 10 kHz
    for f, ref in Y_AITKEN.items():
        got = abs(mg.Y(2j * math.pi * f))
        assert got == pytest.approx(ref, rel=6e-3), (f, got, ref)
    # deep skin: the Mellin asymptote K/sqrt(s) + c_1/s
    K = 6 * L * L * math.sqrt(SIGMA_CU / MU0)
    c1 = -48 * L / (math.pi * MU0)
    s = 2j * math.pi * 1e8
    assert mg.Y(s) == pytest.approx(K / cmath.sqrt(s) + c1 / s, rel=5e-3)


def test_projected_beats_additive_in_the_wall_band(cube):
    """The additive Y_mixed is tens of percent off where the projection is 0.4 %."""
    mesh, mg = cube
    from ngsolve import TaskManager
    with TaskManager():
        lam, tau, g_n, V = A.bulk_foster_via_eigen(mesh, SIGMA_CU, MU0, n_eigen=12)
    K = A.K_SIBC_total(6 * L * L, SIGMA_CU, MU0)
    c1 = -48 * L / (math.pi * MU0)
    s = 2j * math.pi * 1e4
    add = abs(A.Y_mixed(s, lam, tau, g_n, K, c1))
    assert abs(add / Y_AITKEN[1e4] - 1) > 0.2
    assert abs(abs(mg.Y(s)) / Y_AITKEN[1e4] - 1) < 6e-3


def test_rejects_left_half_plane():
    with pytest.raises(ValueError):
        SC.box_layer_1d  # closed forms exist for any t
        SC.BoxMixedGalerkin._t(type("S", (), {"mu": MU0, "sigma": SIGMA_CU})(), -1.0)


def test_multiport_symmetry_and_dc_gram(cube):
    from ngsolve import TaskManager, x, y
    mesh, _ = cube
    with TaskManager():
        mg2 = SC.BoxMixedGalerkin(mesh, SIGMA_CU, MU0, drive_cfs=[1.0, x - L / 2],
                                  n_eigen=12, rule={"ratio": 2.5})
    Y0 = mg2.Y(0)
    assert Y0[0, 0] == pytest.approx(SIGMA_CU * L ** 3, rel=1e-12)
    assert Y0[1, 1] == pytest.approx(SIGMA_CU * L ** 5 / 12, rel=1e-6)
    assert abs(Y0[0, 1]) < 1e-9 * Y0[0, 0]
    for f in (1e3, 1e6):
        Y = mg2.Y(2j * math.pi * f)
        assert Y.shape == (2, 2)
        assert Y[0, 1] == pytest.approx(Y[1, 0], abs=1e-12 * abs(Y[0, 0]))
        assert abs(Y[0, 1]) < 1e-3 * abs(Y[0, 0])       # odd port decouples by symmetry
        assert Y[1, 1].real > 0.0
    a = mg2.alpha(2j * math.pi * 1e6)
    assert a.shape == (2, 2)
    assert 0.9 < (a[0, 0] / mg2.V).real < 1.0


def test_validation_artifact_records_projected_and_time_domain_cases():
    artifact = ROOT / "validation_test/mixed_galerkin/results/mixed_galerkin_results.json"
    data = json.loads(artifact.read_text(encoding="utf-8"))
    time_case = data["cases"]["cube3d_time_domain"]
    edge_case = data["cases"]["box_edge_corner_dofs"]

    assert time_case["script"] == "time_domain/01_cube_step_response.py"
    assert time_case["all_poles_real"]
    assert time_case["residues_nonnegative"]
    assert edge_case["script"] == "cube3d/10_edge_corner_dofs.py"
    assert edge_case["cube"]
