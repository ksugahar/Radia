"""Independent arithmetic and zero-contrast checks for the dual FEM diagnostic."""

from math import isclose

import ngsolve as ng
from netgen.csg import unit_cube
from pytest import raises

from radia.dual_fem_delta import ellipsoid_delta_benchmark, ellipsoid_delta_series, evaluate_delta
from radia.kelvin_material import MU_0


def test_delta_matches_integrated_constant_fields_and_handles_zero_gap():
    with ng.TaskManager():
        mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.6))
        zero = ng.CoefficientFunction((0, 0, 0))
        unit = ng.CoefficientFunction((1, 0, 0))
        regions = {
            "iron": mesh.Materials(".*"),
            "air": mesh.Materials("no_air"),
            "kelvin": mesh.Materials("no_kelvin"),
        }
        case = evaluate_delta(mesh, zero, unit, MU_0, unit, 1,
                              integration_order=8, **regions)
        assert isclose(case["moment_omega"], 1.0)
        assert isclose(case["moment_2_to_1"], 1 / 3)
        assert isclose(case["difference_energy_j"]["iron"], MU_0)
        assert abs(case["delta_j"]) < 1e-14 * MU_0

        vacuum = evaluate_delta(mesh, unit, unit, MU_0, unit, 0,
                                integration_order=8, **regions)
        assert vacuum["moment_gap_j"] == vacuum["delta_j"] == 0
        assert vacuum["relative_delta"] is None


def test_curved_geometry_order_cannot_exceed_field_order():
    with raises(ValueError, match="curve_order"):
        ellipsoid_delta_benchmark(order=2, curve_order=3)


def test_kelvin_geometry_and_midpoint_are_reported():
    with ng.TaskManager():
        result = ellipsoid_delta_series(orders=(2,))
    geometry = result["geometry"]
    assert not geometry["failed_tetrahedra_by_material"]
    assert geometry["min_signed_scaled_jacobian"] > 0
    assert geometry["kelvin_face_pairs"] > 0
    assert geometry["kelvin_curved_match_relative"] <= 1e-10
    row, = result["rows"]
    assert row["semiaxes_m"] == (1.0, 0.8, 0.6)
    assert isclose(row["moment_midpoint"], (row["moment_a"] + row["moment_omega"]) / 2)
    assert set(row["relative_moment_error"]) == {"a", "omega", "2_to_1", "midpoint"}


def test_prolate_requires_affine_mesh_input():
    with raises(ValueError, match="vol_path"):
        ellipsoid_delta_series(shape="prolate")


def test_dual_bracket_reference_freezes_reference_and_applies_r_rule():
    from radia.dual_fem_delta import dual_bracket_reference

    levels = [{"q_lower": 0.90, "q_upper": 1.20}, {"q_lower": 0.999, "q_upper": 1.002},
              {"q_lower": 0.9999, "q_upper": 1.0001}]
    # r = 1.7 at the coarse level: errors -0.1 and +0.17 around 1.0, so 2:1 is closer
    result = dual_bracket_reference(levels, coarse={"q_lower": 0.9, "q_upper": 1.17})
    assert isclose(result["q_ref"], 1.0) and isclose(result["half_width"], 1e-4)
    assert result["converged"] and not result["nesting_violations"]
    assert isclose(result["coarse"]["r"], 1.7)
    assert result["coarse"]["closer"] == "two_to_one"
    # r = 1.2: the midpoint is closer
    assert dual_bracket_reference(levels, coarse={"q_lower": 0.9, "q_upper": 1.12})["coarse"]["closer"] == "midpoint"
    # r = 1.4: the two blends tie, so the comparison is inconclusive
    assert dual_bracket_reference(levels, coarse={"q_lower": 0.9, "q_upper": 1.14})["coarse"]["closer"] == "inconclusive"
    # a finer level outside the coarser bracket is reported and blocks convergence
    bad = dual_bracket_reference([{"q_lower": 0.99, "q_upper": 1.01}, {"q_lower": 0.98, "q_upper": 1.0}])
    assert bad["nesting_violations"] == [1] and not bad["converged"]
    with raises(ValueError):
        dual_bracket_reference([{"q_lower": 1.0, "q_upper": 0.9}, {"q_lower": 0.9, "q_upper": 1.0}])


def test_dual_bracket_reference_rejects_invalid_bounds_before_any_verdict():
    """Guards: reversed or non-finite bounds raise ValueError in both the
    refinement levels and the optional coarse level, before a comparison or
    convergence verdict is formed.  Assumes: equal bounds and no coarse level
    remain valid inputs."""
    from radia.dual_fem_delta import dual_bracket_reference

    nan = float("nan")
    good = [{"q_lower": 9.9, "q_upper": 10.1}, {"q_lower": 9.9, "q_upper": 10.1}]
    bad_levels = [
        [{"q_lower": 10.1, "q_upper": 9.9}, {"q_lower": 9.9, "q_upper": 10.1}],
        [{"q_lower": 9.9, "q_upper": 10.1}, {"q_lower": nan, "q_upper": 10.1}],
        [{"q_lower": 9.9, "q_upper": 10.1}, {"q_lower": 9.9, "q_upper": float("inf")}],
        [{"q_lower": 9.9, "q_upper": 10.1}, {"q_lower": 9.9}],
        [{"q_lower": 9.9, "q_upper": 10.1}, {"q_lower": "x", "q_upper": 10.1}],
    ]
    for levels in bad_levels:
        with raises(ValueError, match=r"levels\[[01]\]"):
            dual_bracket_reference(levels)
    for coarse in ({"q_lower": 12.0, "q_upper": 8.0}, {"q_lower": nan, "q_upper": 10.5},
                   {"q_lower": 9.5, "q_upper": nan}, {"q_upper": 10.5}):
        with raises(ValueError, match="coarse"):
            dual_bracket_reference(good, coarse=coarse)
    # valid edge cases still produce a verdict
    equal = dual_bracket_reference([{"q_lower": 1.0, "q_upper": 1.0}, {"q_lower": 1.0, "q_upper": 1.0}])
    assert equal["half_width"] == 0.0 and "coarse" not in equal
    with_equal_coarse = dual_bracket_reference(good, coarse={"q_lower": 9.5, "q_upper": 9.5})
    assert with_equal_coarse["coarse"]["closer"] in {"midpoint", "two_to_one", "inconclusive"}


def test_dual_bracket_reference_rejects_bad_fraction_and_unrepresentable_values():
    """Guards: a non-finite or negative convergence_fraction raises; finite
    bounds whose derived reference overflows raise instead of returning
    converged=True; an unconverged reference never reports a closer blend and
    keeps the would-be answer as provisional_closer.  Assumes: the default
    fraction 0.1 and the r < 1.4 rule are unchanged."""
    from radia.dual_fem_delta import dual_bracket_reference

    unconverged = [{"q_lower": 0.0, "q_upper": 30.0}, {"q_lower": 9.9, "q_upper": 10.1}]
    coarse = {"q_lower": 9.0, "q_upper": 40.0}
    for fraction in (float("inf"), float("nan"), -0.1, "x", None):
        with raises(ValueError, match="convergence_fraction"):
            dual_bracket_reference(unconverged, coarse=coarse, convergence_fraction=fraction)
    huge = [{"q_lower": 1e308, "q_upper": 1e308}, {"q_lower": 1e308, "q_upper": 1e308}]
    with raises(ValueError, match="not representable"):
        dual_bracket_reference(huge)
    with raises(ValueError, match="not representable"):
        dual_bracket_reference([{"q_lower": 9.9, "q_upper": 10.1}] * 2,
                               coarse={"q_lower": -1e308, "q_upper": 1e308})
    result = dual_bracket_reference(unconverged, coarse=coarse)
    assert result["converged"] is False, result
    assert result["coarse"]["closer"] == "inconclusive", result["coarse"]
    assert result["coarse"]["provisional_closer"] == "two_to_one", result["coarse"]
    assert result["last_step_drift"] == {"lower": 9.9, "upper": -19.9}, result["last_step_drift"]


def test_dual_bracket_reference_rejects_overflowing_and_lossy_inputs():
    """Guards: integers too large for binary64 raise ValueError (not
    OverflowError) in the fraction, the levels and the coarse level; a finite
    convergence budget overflow raises; a reversed pair that rounds to equal
    floats is still rejected.  Assumes: ordinary JSON integers reach the entry
    through the MCP dictionary unchanged."""
    from radia.dual_fem_delta import dual_bracket_reference

    good = [{"q_lower": 9.9, "q_upper": 10.1}, {"q_lower": 9.9, "q_upper": 10.1}]
    huge = 10 ** 400
    with raises(ValueError, match="convergence_fraction"):
        dual_bracket_reference(good, convergence_fraction=huge)
    with raises(ValueError, match=r"levels\[0\]"):
        dual_bracket_reference([{"q_lower": huge, "q_upper": huge}] + good[1:])
    with raises(ValueError, match="coarse"):
        dual_bracket_reference(good, coarse={"q_lower": huge, "q_upper": 8})
    with raises(ValueError, match="convergence_budget"):
        dual_bracket_reference(good, coarse={"q_lower": 0.0, "q_upper": 1e308},
                               convergence_fraction=1e308)
    lossy = {"q_lower": 2 ** 53 + 1, "q_upper": 2 ** 53}
    assert float(lossy["q_lower"]) == float(lossy["q_upper"])
    with raises(ValueError, match=r"levels\[0\].*q_lower <= q_upper"):
        dual_bracket_reference([lossy, lossy])
    with raises(ValueError, match="coarse.*q_lower <= q_upper"):
        dual_bracket_reference(good, coarse=lossy)
    ordered = dual_bracket_reference([{"q_lower": 2 ** 53, "q_upper": 2 ** 53 + 1}] * 2)
    assert ordered["half_width"] == 0.0, ordered
