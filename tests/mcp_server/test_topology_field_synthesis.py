"""Golden lock for the topology_optimization `field_synthesis` topic.

Two jobs:
  1. SYMBOLICALLY re-verify the analytic PM-multipole inverse formula so the
     knowledge text can never silently rot (Hidaka single-mode + Sugahara
     general radially-varying form).  Both must give Laplacian(A_z) - rotM = 0.
  2. Check the topic is wired into get_applications_documentation and carries
     the load-bearing facts (the 1/(n^2-1) coefficient, the n=1 degeneracy,
     and the pointer to the streamfunction coil branch).
"""

import pytest

from radia_mcp.topology_optimization.applications_knowledge import (
    get_applications_documentation,
)


def test_pm_multipole_inverse_is_exact():
    """A_z = mu0 M_rn r cos(n(th+th0))/(n^2-1) solves the polar Poisson eqn."""
    sp = pytest.importorskip("sympy")
    r, th, th0, Mrn = sp.symbols("r theta theta_0 Mrn", positive=True)
    n = sp.symbols("n", integer=True, positive=True)

    # Hidaka single-mode, r-independent amplitude (mu0 = 1, lab convention)
    Mr = Mrn * (sp.cos(n * th0) * sp.sin(n * th)
                + sp.sin(n * th0) * sp.cos(n * th)) / n
    rotMr = -sp.diff(Mr, th) / r
    Az = -(r * Mrn) * (-sp.cos(n * th0) * sp.cos(n * th)
                       + sp.sin(n * th0) * sp.sin(n * th)) / (n ** 2 - 1)
    lap = sp.diff(r * sp.diff(Az, r), r) / r + sp.diff(Az, th, 2) / r ** 2
    assert sp.simplify(lap - rotMr) == 0


def test_pm_multipole_general_radial_is_exact():
    """Sugahara general radially-varying form (M_r ~ r^{+/-n}) is also exact."""
    sp = pytest.importorskip("sympy")
    r, th = sp.symbols("r theta", positive=True)
    A, B, C, D = sp.symbols("A B C D", real=True)
    n = sp.symbols("n", integer=True, positive=True)

    Mr = (A * r ** n * sp.cos(n * th) + B * r ** (-n) * sp.cos(n * th)
          + C * r ** n * sp.sin(n * th) + D * r ** (-n) * sp.sin(n * th))
    rotMr = -sp.diff(Mr, th) / r
    Az = (-(n / (2 * n + 1)) * C * r ** (n + 1) * sp.cos(n * th)
          + (n / (2 * n - 1)) * D * r ** (-n + 1) * sp.cos(n * th)
          + (n / (2 * n + 1)) * A * r ** (n + 1) * sp.sin(n * th)
          + (n / (-2 * n + 1)) * B * r ** (-n + 1) * sp.sin(n * th))
    lap = sp.diff(r * sp.diff(Az, r), r) / r + sp.diff(Az, th, 2) / r ** 2
    assert sp.simplify(lap - rotMr) == 0


def test_field_synthesis_topic_wired_and_loadbearing():
    doc = get_applications_documentation("field_synthesis")
    assert doc.strip()
    # the modal coefficient and the degeneracy must be stated
    assert "n^2 - 1" in doc
    assert "n = 1" in doc and "degenerac" in doc.lower()
    # must point to the shipped coil branch, not reimplement it
    assert "streamfunction" in doc.lower()
    # "all" must concatenate motor + field synthesis
    all_doc = get_applications_documentation("all")
    assert "IPM motor" in all_doc and "Field synthesis" in all_doc


def test_unknown_topic_message():
    out = get_applications_documentation("nonsense")
    assert "Unknown topic" in out and "field_synthesis" in out


def test_outer_loop_topic_wired():
    doc = get_applications_documentation("outer_loop")
    assert doc.strip()
    assert "Nelder-Mead" in doc and "fminsearchbnd" in doc
    # "all" must include it alongside motor + field synthesis
    all_doc = get_applications_documentation("all")
    assert "Outer-loop optimizers" in all_doc


def test_nelder_mead_and_bound_transform():
    """Ground the documented numbers: NM->(1,1); fminsearchbnd LB=2 transform->(2,4)."""
    np = pytest.importorskip("numpy")
    opt = pytest.importorskip("scipy.optimize")
    rosen = lambda v: (1 - v[0]) ** 2 + 105 * (v[1] - v[0] ** 2) ** 2
    r = opt.minimize(rosen, [3.0, 3.0], method="Nelder-Mead",
                     options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 5000})
    assert np.allclose(r.x, [1.0, 1.0], atol=1e-3)
    # D'Errico one-sided (quadratic) transform: x = LB + t^2, LB=(2,2)
    LB = np.array([2.0, 2.0])
    obj_t = lambda t: rosen(LB + t ** 2)
    rt = opt.minimize(obj_t, np.sqrt([1.0, 1.0]), method="Nelder-Mead",
                      options={"xatol": 1e-8, "fatol": 1e-10, "maxiter": 5000})
    x = LB + rt.x ** 2
    assert np.allclose(x, [2.0, 4.0], atol=1e-2) and abs(rt.fun - 1.0) < 1e-2
