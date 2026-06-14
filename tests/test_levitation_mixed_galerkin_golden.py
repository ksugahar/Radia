"""Golden tests for the radia-levitation Mixed-Galerkin polarizability stack.

Locks the numerical invariants of the alpha(s) pipeline (Phase L-3):

  * the closed-form wedge function W(alpha) = (4/pi) cot(alpha/2),
  * the planar SIBC anchor K_SIBC = S sqrt(sigma/mu),
  * the cuboid edge correction c_1 = -16(a+b+c)/(pi mu) (CAD-direct),
  * the perfect-conductor high-frequency limit alpha(s)->V,
  * the Simulink/MATLAB state-space LTI vs the closed-form Y(s),
  * an end-to-end NGSolve cube alpha(s) sweep.

Pure-algebra tests always run.  Tests needing netgen.occ / ngsolve
`importorskip` out when those are absent (e.g. a minimal-dep CI box).

Scope note: the Mixed-Galerkin Y(s) carries a high-frequency Mellin tail
K_SIBC/sqrt(s) + c_1/s which DIVERGES as s->0, so alpha(s)/V is only the
physical polarizability in the mid-to-high band, NOT at DC.  The golden
limit locked here is therefore the PEC limit alpha(s->inf)->V, not the
DC limit.
"""
from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

from radia.levitation.mixed_galerkin import alpha as A
from radia.levitation.simulink import export as EX

MU0 = 4.0 * math.pi * 1e-7
SIGMA_CU = 5.8e7


# ---------------------------------------------------------------------------
# Pure-algebra invariants (no external solver)
# ---------------------------------------------------------------------------

def test_wedge_function_special_values():
    """W(pi/2) = 4/pi (right angle), W(pi) = 0 (flat), W(3pi/2) = -4/pi (re-entrant)."""
    assert A.wedge_function(math.pi / 2) == pytest.approx(4.0 / math.pi, rel=1e-12)
    assert A.wedge_function(math.pi * (1 - 1e-9)) == pytest.approx(0.0, abs=1e-5)
    assert A.wedge_function(3 * math.pi / 2) == pytest.approx(-4.0 / math.pi, rel=1e-12)


def test_K_SIBC_total_formula():
    """K_SIBC = S sqrt(sigma/mu)."""
    S = 6.0e-4
    assert A.K_SIBC_total(S, SIGMA_CU, MU0) == pytest.approx(
        S * math.sqrt(SIGMA_CU / MU0), rel=1e-12
    )


def test_c1_polyhedral_cube_closed_form():
    """c_1 = -(1/mu) sum L_e W(alpha_e); a cube has 12 edges of length L at pi/2.

    Closed form: c_1 = -(1/mu) * 12 L * (4/pi) = -48 L / (pi mu)
                     = -16 (a+b+c) / (pi mu)  with a=b=c=L (3L total / axis).
    """
    L = 5e-3
    edges = [(L, math.pi / 2)] * 12
    c1 = A.c1_polyhedral(edges, MU0)
    closed = -16.0 * (3 * L) / (math.pi * MU0)
    assert c1 == pytest.approx(closed, rel=1e-12)


def test_alpha_pec_high_frequency_limit():
    """alpha(s)/V -> 1 as s -> inf (perfect-conductor flux exclusion)."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0])
    tau = MU0 * SIGMA_CU / lam
    g_n = np.array([SIGMA_CU * V * 0.30, SIGMA_CU * V * 0.10, SIGMA_CU * V * 0.05])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    c1 = -6.08e4
    a_hi = A.alpha_from_Y(A.Y_mixed(1j * 2 * math.pi * 1e12, lam, tau, g_n, K, c1),
                          V, SIGMA_CU)
    assert (a_hi / V).real == pytest.approx(1.0, abs=2e-3)
    assert abs((a_hi / V).imag) < 1e-2


def test_state_space_lti_matches_closed_form_Y():
    """ss(A,B,C,D) transfer function reproduces Y_mixed(s) to <1% over 1 Hz-10 kHz."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0])
    tau = MU0 * SIGMA_CU / lam
    g_n = np.array([SIGMA_CU * V * 0.30, SIGMA_CU * V * 0.10, SIGMA_CU * V * 0.05])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    c1 = -6.08e4

    A_, B_, C_, D_, n_f, n_w, n_i = EX.build_state_space(
        g_n, tau, V, SIGMA_CU, K, c1, n_warburg_rungs=40
    )
    assert n_f == len(g_n)
    assert n_i == 1

    eye = np.eye(A_.shape[0])
    for f in (1.0, 1e2, 1e4):
        s = 1j * 2 * math.pi * f
        y_ss = (C_ @ np.linalg.solve(s * eye - A_, B_) + D_)[0, 0]
        y_cf = A.Y_mixed(s, lam, tau, g_n, K, c1)
        assert abs(y_ss - y_cf) / abs(y_cf) < 1e-2


def test_diffusive_quadrature_residues_positive():
    """Warburg diffusive-quadrature rungs are positive and log-spaced."""
    K = 6.8e6
    xi, r = EX.diffusive_quadrature(K, n_aux=30)
    assert len(xi) == len(r) == 30
    assert np.all(r > 0)
    assert np.all(np.diff(np.log(xi)) > 0)  # strictly log-increasing poles


# ---------------------------------------------------------------------------
# CAD-direct edge extraction (needs netgen.occ)
# ---------------------------------------------------------------------------

def test_cube_cad_topology_c1():
    """CAD-direct: a Box has exactly 12 edges and c_1 == -16(3L)/(pi mu)."""
    occ = pytest.importorskip("netgen.occ")
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    c1, L_total, n_edges = CE.cad_topology_c1(box, MU0)

    assert n_edges == 12
    assert L_total == pytest.approx(12 * L, rel=1e-9)
    closed = -16.0 * (3 * L) / (math.pi * MU0)
    assert c1 == pytest.approx(closed, rel=1e-6)


# ---------------------------------------------------------------------------
# End-to-end NGSolve cube alpha(s) (needs ngsolve; fast, ~1 s)
# ---------------------------------------------------------------------------

def test_cube_alpha_sweep_end_to_end():
    """Coarse Cu cube: V exact, PEC high-f limit, c_1 closed form.

    Geometry-exact quantities (V, c_1) are locked tightly; the bulk-Foster
    PEC limit alpha(1 GHz)/V is locked to a band (mesh/eigen-count
    dependent, approaches 1.0 from below).
    """
    pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")
    from ngsolve import Mesh, TaskManager
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    for f in box.faces:
        f.name = "outer"
    ng = occ.OCCGeometry(box).GenerateMesh(maxh=L / 4)
    mesh = Mesh(ng)

    with TaskManager():
        lam, tau, g_n, V = A.bulk_foster_via_eigen(mesh, SIGMA_CU, MU0, n_eigen=40)

    # Geometry-exact invariants
    assert V == pytest.approx(L**3, rel=2e-2)          # 125 mm^3
    c1, _, n_edges = CE.cad_topology_c1(box, MU0)
    assert n_edges == 12
    assert c1 == pytest.approx(-16.0 * (3 * L) / (math.pi * MU0), rel=1e-6)

    # Foster spectrum sane
    assert len(lam) == 40
    assert np.all(lam > 0)
    assert np.all(g_n >= 0)

    # PEC high-frequency limit alpha(1 GHz)/V in [0.95, 1.02]
    K = A.K_SIBC_total(6 * L * L, SIGMA_CU, MU0)
    s = 1j * 2 * math.pi * 1e9
    a = A.alpha_from_Y(A.Y_mixed(s, lam, tau, g_n, K, c1), V, SIGMA_CU)
    assert 0.95 < (a / V).real < 1.02
