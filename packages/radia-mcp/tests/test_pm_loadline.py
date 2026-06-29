r"""PM operating point / load line in a gapped iron circuit -- #42.

A magnet in a window frame drives flux across an air gap; the load-line gap field
B_gap = Br l_m/(l_m + mu_rec g + l_fe/mu_r) is the leakage-free UPPER bound. Tool-independent
gate: the FE B_gap is BELOW the load line (window leakage), the deficit GROWS with the gap,
and B_gap falls monotonically with the gap (more demagnetisation). A radia self-regression
pins the values (an independent 2D FE solver matches them to ~0.3 %, recorded internally).
The PM counterpart of the coil-driven #27 (force) / #39 (inductance)."""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.solve import (
    pm_circuit_loadline_gap_field,
    pm_circuit_loadline_operating_point,
    pm_loadline_demag_risk_summary,
    pm_temperature_demag_sweep_summary,
    MU0,
)

W, ww, b = 0.24, 0.12, 0.06
MU_R, RBOX = 2000.0, 0.60
BR, MU_REC, LM = 1.2, 1.0, 0.06

# radia self-regression references (tesla) -- an independent 2D FE solver matches to ~0.3%
_REF = {0.002: 0.9914, 0.004: 0.8700, 0.008: 0.7109}


def test_loadline_formula():
    # B_gap = Br l_m/(l_m + mu_rec g + l_fe/mu_r); larger gap -> lower B; mu_r->inf drops iron term
    g, lfe = 0.004, 0.6
    B = pm_circuit_loadline_gap_field(BR, LM, g, lfe, MU_R, MU_REC)
    assert math.isclose(B, BR*LM/(LM + MU_REC*g + lfe/MU_R), rel_tol=1e-12)
    assert pm_circuit_loadline_gap_field(BR, LM, 0.008, lfe, MU_R) < B      # bigger gap -> lower
    # zero gap, ideal iron -> B_gap -> Br (magnet fully short-circuited)
    assert math.isclose(pm_circuit_loadline_gap_field(BR, LM, 0.0, 0.0, 1e12), BR, rel_tol=1e-9)


def test_loadline_operating_point_identities():
    g, lfe, mu_rec, hknee = 0.004, 0.6, 1.05, -6.0e5
    op = pm_circuit_loadline_operating_point(BR, LM, g, lfe, MU_R, mu_rec, H_knee=hknee)
    expected_b = pm_circuit_loadline_gap_field(BR, LM, g, lfe, MU_R, mu_rec)
    expected_pc = LM / (g + lfe / (mu_rec * MU_R))

    assert op["B_gap_T"] == pytest.approx(expected_b)
    assert op["H_m_A_per_m"] == pytest.approx((expected_b - BR) / (MU0 * mu_rec))
    assert op["permeance_coefficient"] == pytest.approx(expected_pc)
    assert op["B_from_permeance_coefficient_T"] == pytest.approx(expected_b)
    assert op["B_identity_abs_error_T"] < 1.0e-15
    assert op["demag_margin_A_per_m"] == pytest.approx(op["H_m_A_per_m"] - hknee)
    assert op["safe_against_knee"] is True


def test_loadline_operating_point_gap_sweep_monotone():
    gaps = (0.001, 0.002, 0.004, 0.008)
    ops = [pm_circuit_loadline_operating_point(BR, LM, g, 0.6, MU_R, MU_REC) for g in gaps]
    bs = [op["B_gap_T"] for op in ops]
    hs = [op["H_m_A_per_m"] for op in ops]
    pcs = [op["permeance_coefficient"] for op in ops]

    assert all(a > b for a, b in zip(bs, bs[1:]))
    assert all(a > b for a, b in zip(hs, hs[1:]))  # H_m becomes more negative.
    assert all(a > b for a, b in zip(pcs, pcs[1:]))


def test_continuous_loop_pm_loadline_demag_margin_gate():
    params = {
        "Br": 1.2,
        "magnet_len": 0.004,
        "iron_path": 0.08,
        "mu_r": 1000.0,
        "mu_rec": 1.05,
        "H_knee": -4.5e5,
    }
    rows = [
        pm_circuit_loadline_operating_point(gap=gap, **params)
        for gap in (0.0005, 0.001, 0.002, 0.004, 0.008)
    ]

    assert rows[0]["B_gap_T"] == pytest.approx(1.042345276872964)
    assert rows[1]["permeance_coefficient"] == pytest.approx(3.7168141592920354)
    assert rows[2]["H_m_A_per_m"] == pytest.approx(-320811.6282388644)
    assert rows[3]["demag_margin_A_per_m"] == pytest.approx(-20105.698021609453)
    assert rows[-1]["B_identity_abs_error_T"] < 1.0e-14
    assert [row["safe_against_knee"] for row in rows] == [True, True, True, False, False]


def test_pm_loadline_demag_risk_summary_marks_safe_prefix():
    params = {
        "Br": 1.068,
        "magnet_len": 0.004,
        "iron_path": 0.08,
        "mu_r": 1000.0,
        "mu_rec": 1.05,
        "H_knee": -4.95e5,
    }
    rows = [
        {"gap_m": gap, **pm_circuit_loadline_operating_point(gap=gap, **params)}
        for gap in (0.0005, 0.001, 0.002, 0.004, 0.008)
    ]

    summary = pm_loadline_demag_risk_summary(rows, axis_key="gap_m")

    assert summary["schema"] == "radia-ngsolve.pm-loadline-demag-risk-summary.v1"
    assert summary["risk_label"] == "red"
    assert summary["safe_prefix_count"] == 4
    assert summary["first_unsafe_index"] == 4
    assert summary["first_unsafe_gap_m"] == pytest.approx(0.008)
    assert summary["largest_safe_gap_m"] == pytest.approx(0.004)
    assert summary["minimum_demag_margin_A_per_m"] == pytest.approx(-54988.18063331157)
    assert summary["checks"]["safe_region_is_prefix"] is True
    assert summary["checks"]["margin_monotone_decreasing"] is True


def test_continuous_loop_pm_temperature_demag_sweep_gate():
    sweep = pm_temperature_demag_sweep_summary(
        Br_20C=1.2,
        H_knee_20C=-9.0e5,
        temperature_C=120.0,
        magnet_len=0.004,
        gaps=(0.0005, 0.001, 0.002, 0.004, 0.008),
        iron_path=0.08,
        mu_r=1000.0,
        mu_rec=1.05,
    )

    assert sweep["schema"] == "radia-ngsolve.pm-temperature-demag-sweep.v1"
    assert sweep["Br_hot_T"] == pytest.approx(1.068)
    assert sweep["H_knee_hot_A_per_m"] == pytest.approx(-495000.0)
    assert sweep["first_unsafe_gap_m"] == pytest.approx(0.008)
    assert sweep["safe_prefix_count"] == 4
    assert sweep["risk_label"] == "red"
    assert sweep["risk_summary"]["largest_safe_gap_m"] == pytest.approx(0.004)
    assert all(sweep["checks"].values())
    assert [row["safe_against_knee"] for row in sweep["rows"]] == [True, True, True, True, False]
    assert sweep["rows"][0]["demag_margin_A_per_m"] == pytest.approx(388659.7126007991)
    assert sweep["rows"][-1]["demag_margin_A_per_m"] == pytest.approx(-54988.18063331157)


def test_continuous_loop_slot46_elf_hbrm_hbcn_loadline_contract():
    sweep = pm_temperature_demag_sweep_summary(
        Br_20C=1.2,
        H_knee_20C=-9.0e5,
        temperature_C=120.0,
        magnet_len=0.004,
        gaps=(0.0005, 0.001, 0.002, 0.004, 0.008),
        iron_path=0.08,
        mu_r=1000.0,
        mu_rec=1.05,
    )

    assert sweep["Br_hot_T"] == pytest.approx(1.068)
    assert sweep["H_knee_hot_A_per_m"] == pytest.approx(-495000.0)
    assert sweep["safe_prefix_count"] == 4
    assert sweep["first_unsafe_gap_m"] == pytest.approx(0.008)
    assert sweep["risk_summary"]["largest_safe_gap_m"] == pytest.approx(0.004)
    assert sweep["minimum_demag_margin_A_per_m"] == pytest.approx(-54988.18063331157)
    assert sweep["risk_label"] == "red"
    assert all(sweep["checks"].values())

    expected = [
        (0.0005, 0.9276872964169381, 388659.7126007991, True),
        (0.001, 0.8327485380116959, 316707.4616903156, True),
        (0.002, 0.6912621359223301, 209477.6508674107, True),
        (0.004, 0.5159420289855072, 76605.92876076762, True),
        (0.008, 0.3423076923076923, -54988.18063331157, False),
    ]
    for row, (gap, b_gap, margin, safe) in zip(sweep["rows"], expected):
        assert row["gap_m"] == pytest.approx(gap)
        assert row["B_gap_T"] == pytest.approx(b_gap)
        assert row["demag_margin_A_per_m"] == pytest.approx(margin)
        assert row["safe_against_knee"] is safe


def _rect(wp, x0, y0, w, h):
    return wp.MoveTo(x0, y0).Rectangle(w, h).Face()


def _Bgap_fe(g):
    from ngsolve import Mesh, CoefficientFunction, grad, Integrate, dx, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane, Glue, X, Y
    from radia_mcp.radia_ngsolve.solve import solve_planar_magnetostatic, NU0
    big = WorkPlane().MoveTo(-W/2, -W/2).Rectangle(W, W).Face()
    window = WorkPlane().MoveTo(-ww/2, -ww/2).Rectangle(ww, ww).Face()
    gap = WorkPlane().MoveTo(-g/2, ww/2).Rectangle(g, b).Face()
    mag = WorkPlane().MoveTo(-W/2, -LM/2).Rectangle((W-ww)/2, LM).Face()
    iron = big - window - gap - mag
    gap.faces.name = "gap"; mag.faces.name = "magnet"; iron.faces.name = "iron"
    box = WorkPlane().MoveTo(-RBOX, -RBOX).Rectangle(2*RBOX, 2*RBOX).Face()
    air = box - iron - gap - mag; air.faces.name = "air"
    air.edges.Max(X).name = "outer"; air.edges.Min(X).name = "outer"
    air.edges.Max(Y).name = "outer"; air.edges.Min(Y).name = "outer"
    mesh = Mesh(OCCGeometry(Glue([air, iron, gap, mag]), dim=2).GenerateMesh(maxh=0.04))
    nu = mesh.MaterialCF({"iron": NU0/MU_R, "magnet": NU0/MU_REC}, default=NU0)
    with TaskManager():
        gfu = solve_planar_magnetostatic(mesh, nu, magnets={"magnet": (BR/MU0, 90.0)},
                                         order=3, dirichlet="outer")
    B = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    Agap = Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*dx, mesh)
    return abs(Integrate(mesh.MaterialCF({"gap": 1.0}, default=0.0)*B[0]*dx, mesh) / Agap)
