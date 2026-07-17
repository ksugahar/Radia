"""Fast contracts for the production C++ planar HDiv-VIM image operator."""

import numpy as np
import ngsolve as ng
import pytest
from ngsolve.meshes import MakeStructured2DMesh

from radia import vim


def _full_square():
    return MakeStructured2DMesh(
        quads=True, nx=2, ny=2, mapping=lambda x, y: (2.0*x - 1.0, 2.0*y - 1.0))


def _right_half_square():
    return MakeStructured2DMesh(
        quads=True, nx=1, ny=2, mapping=lambda x, y: (x, 2.0*y - 1.0))


@pytest.mark.parametrize(
    ("image", "field"),
    [("+x", (0.0, 1.0)), ("-x", (1.0, 0.0))],
)
@pytest.mark.parametrize("order", [1, 2])
def test_planar_quad_image_solve_and_field_match_full_to_roundoff(
        image, field, order):
    """Both reflection parities reproduce the matching explicit full mesh."""
    applied = ng.CoefficientFunction(field)
    probes = np.array([[2.0, 0.0], [1.5, 0.3], [0.2, 1.5]])
    with ng.TaskManager():
        full = vim.Solve(
            _full_square(), order=order, mu_r=100.0, H_ext=applied, tol=1e-13)
        half = vim.Solve(
            _right_half_square(), order=order, mu_r=100.0,
            H_ext=applied, image=image, tol=1e-13)
        full_field = full["body"].H_at(probes, full["m"])
        half_field = half["body"].H_at(probes, half["m"])

    relative = np.linalg.norm(full_field - half_field) / np.linalg.norm(full_field)
    # BDM2 adds two local reference-basis transformations before the C++ solve;
    # its 12-20 eps full/reduced spread is still roundoff, not a discretization
    # tolerance.  BDM1 keeps the original strict 10 eps gate.
    roundoff_factor = 10.0 if order == 1 else 32.0
    assert relative < roundoff_factor * np.finfo(float).eps
    expected_constraints = 2*(order+1) if image == "+x" else 0
    assert half["body"].G.constraint_count == expected_constraints
    assert half["field_evaluator_stats"]["image_count"] == 1


def test_charge_gram_exposes_only_persistent_operator_bindings():
    """The topology-resending pybind surface is intentionally not backward compatible."""
    with ng.TaskManager():
        result = vim.Solve(
            _right_half_square(), mu_r=100.0,
            H_ext=ng.CoefficientFunction((1.0, 0.0)), image="-x")
    gram = result["body"].G
    assert gram.operator_configured
    for removed in (
        "solve_linear_material_auto_prec",
        "solve_linear_material_mass_riesz",
        "apply_demag",
        "apply_mass_riesz",
    ):
        assert not hasattr(gram, removed)
