"""Golden: 2D LINEAR hodograph feasibility of a bending magnet.

Locks docs/clebsch_hodograph/hodograph_feasibility_2d.{ipynb,py}:

  (1) feasibility law  d* = (2/pi) h  (edge no sharper than ~0.64 x gap): the demanded
      field is realizable by an iron pole at gap h iff the tanh continuation's singularity
      pi*d/2 clears the gap;
  (2) inverse design + independent ngsolve linear-FEM verification (manufactured solution):
      the read-off equipotential pole reproduces the demanded mid-plane field to ~0.01 %,
      the FEM phi matches the analytic phi, and the FEM equipotential matches the read-off
      pole to ~1e-6.
"""
import math
import sys
from pathlib import Path

import pytest

DOCDIR = Path(__file__).resolve().parents[2] / "docs" / "clebsch_hodograph"
sys.path.insert(0, str(DOCDIR))


def test_feasibility_law_analytic():
    import hodograph_feasibility_2d as hf
    # d* = (2/pi) h exactly
    assert abs(hf.d_star(1.0) - 2.0 / math.pi) < 1e-12
    assert abs(hf.d_star(2.0) - 4.0 / math.pi) < 1e-12
    # classification: singularity pi*d/2 vs gap
    tab = {round(r["d"], 3): r["feasible"] for r in
           hf.feasibility_table([1.4, 1.0, 0.8, 0.5, 0.35])}
    assert tab[1.4] and tab[1.0] and tab[0.8]              # d > d* ~ 0.637
    assert not tab[0.5] and not tab[0.35]                  # d < d*
    # the demo edge d=0.9 is feasible
    assert hf.y_sing(hf.D) > hf.H


def test_vonmises_partial_hodograph():
    import hodograph_feasibility_2d as hf
    # partial von-Mises (x, s=-phi) chart: the free-boundary iron pole becomes a fixed top edge
    vm = hf.vonmises_chart(nX=161, nS=121)
    # (a) single-valued diffeomorphism, no fold: Jacobian y_s > 0 everywhere
    assert vm["single_valued"] and vm["jac_min"] > 0.3, vm
    # (b) the closed-form map satisfies the von-Mises PDE (residual is FD-limited)
    assert vm["pde_resid"] < 3e-3, vm
    # (c) residual drops ~O(h^2) under refinement -> the PDE is right, not a coincidence
    vm_c = hf.vonmises_chart(nX=81, nS=61)
    assert vm["pde_resid"] < 0.6 * vm_c["pde_resid"], (vm["pde_resid"], vm_c["pde_resid"])
    # (d) the chart's top edge s=s0 IS the equipotential iron pole
    assert vm["pole_match"] < 1e-6, vm


@pytest.mark.slow
def test_fem_verifies_inverse_design():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import hodograph_feasibility_2d as hf

    res = hf.fem_verify(order=4)
    assert res["feasible"] is True, res
    # (a) FEM reproduces the demanded mid-plane field
    assert res["By_rel_err"] < 1e-3, res                  # ~1e-4 (0.01 %)
    # (b) FEM phi == analytic phi (independent discretization)
    assert res["phi_L2_err"] < 1e-5, res                  # ~1e-8
    # (c) FEM equipotential == read-off pole; pole passes through (0, h)
    assert res["pole_match_err"] < 5e-4, res              # ~1e-6
    assert abs(res["pole_gap_at_x0"] - hf.H) < 5e-3, res
