"""Native Kelvin source contracts for compact Radia currents."""

from __future__ import annotations

import math

import numpy as np
import pytest


MU_0 = 4.0e-7 * math.pi


@pytest.mark.parametrize(
    "factory_name",
    (
        "KelvinRadiaVectorPotential",
        "KelvinRadiaFluxDensity",
        "KelvinRadiaFieldStrength",
        "KelvinRadiaScalarPotential",
    ),
)
def test_native_kelvin_sources_reject_nonfinite_centres(factory_name):
    import radia as rad

    factory = getattr(rad, factory_name)
    with pytest.raises(ValueError, match="kelvin_center must be finite"):
        factory(1, (math.nan, 0.0, 0.0), 1.0, (0.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="physical_center must be finite"):
        factory(1, (0.0, 0.0, 0.0), 1.0, (0.0, math.inf, 0.0))


def _pullback(values, point, center, radius, *, twisted):
    """Evaluate the analytic Kelvin transform at one computational point."""
    delta = np.asarray(point, dtype=float) - np.asarray(center, dtype=float)
    rho2 = float(np.dot(delta, delta))
    scale = radius * radius / rho2
    normal = delta / math.sqrt(rho2)
    reflected = values - 2.0 * np.dot(values, normal) * normal
    return (-scale if twisted else scale) * reflected, scale


def test_native_kelvin_h_source_obeys_hodge_contract_and_batch_evaluation():
    """The C++ H source has the twisted 1-form sign and matches B = mu H."""
    import ngsolve as ng
    from netgen.csg import unit_cube
    import radia as rad

    kelvin_center = np.array([0.5, 0.5, 0.5])
    physical_center = np.zeros(3)
    radius = 0.2
    points = np.array([
        [0.62, 0.53, 0.49],
        [0.37, 0.54, 0.52],
        [0.53, 0.36, 0.48],
    ])
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.5))

    rad.UtiDelAll()
    try:
        coil = rad.ObjRecCur([0.0, 0.0, 0.0], [0.02, 0.20, 0.02],
                              [0.0, 2.0e6, 0.0])
        h_kelvin = rad.KelvinRadiaFieldStrength(
            coil, tuple(kelvin_center), radius, tuple(physical_center))
        b_kelvin = rad.KelvinRadiaFluxDensity(
            coil, tuple(kelvin_center), radius, tuple(physical_center))

        for point in points:
            delta = point - kelvin_center
            rho2 = float(np.dot(delta, delta))
            mapped = physical_center + radius * radius / rho2 * delta
            h_physical = np.asarray(rad.Fld(coil, "h", [mapped.tolist()])[0])
            b_physical = np.asarray(rad.Fld(coil, "b", [mapped.tolist()])[0])
            expected_h, scale = _pullback(
                h_physical, point, kelvin_center, radius, twisted=True)
            expected_b, _ = _pullback(
                b_physical, point, kelvin_center, radius, twisted=True)
            # A 2-form receives one additional conformal factor.
            expected_b *= scale

            got_h = np.asarray(h_kelvin(mesh(*point)), dtype=float)
            got_b = np.asarray(b_kelvin(mesh(*point)), dtype=float)
            # Radia's scalar field API and its native batch API differ in
            # low-order rounding of near-zero components.  Assess the vector
            # transform against the physical field norm, not a tiny component.
            h_relative_error = np.linalg.norm(got_h - expected_h) / np.linalg.norm(expected_h)
            b_relative_error = np.linalg.norm(got_b - expected_b) / np.linalg.norm(expected_b)
            assert h_relative_error < 1.0e-10
            assert b_relative_error < 1.0e-10
            assert np.linalg.norm(got_h) > 1.0e-12
            assert np.allclose(got_b, MU_0 * scale * got_h,
                               rtol=2e-11, atol=2e-18)

        # GridFunction.Set invokes the coefficient's integration-rule overload
        # and therefore locks the native batch-evaluation path too.
        with ng.TaskManager():
            source_fes = ng.VectorL2(mesh, order=2)
            source = ng.GridFunction(source_fes)
            source.Set(h_kelvin)
        assert np.isfinite(source.vec.FV().NumPy()).all()
    finally:
        rad.UtiDelAll()


def test_native_kelvin_scalar_potential_has_the_matching_twisted_zero_form():
    """A valid Radia scalar trace keeps H' = -grad(phi') after Kelvin inversion."""
    import ngsolve as ng
    from netgen.csg import unit_cube
    import radia as rad

    kelvin_center = np.array([0.5, 0.5, 0.5])
    physical_center = np.zeros(3)
    radius = 0.2
    point = np.array([0.64, 0.54, 0.47])
    mesh = ng.Mesh(unit_cube.GenerateMesh(maxh=0.5))

    rad.UtiDelAll()
    try:
        # A full annular arc is one of the native Radia source types whose
        # scalar potential is locally single-valued away from its current body.
        coil = rad.ObjArcCur(
            [0.0, 0.0, 0.0], [0.08, 0.10], [0.0, 2.0 * math.pi], 0.02,
            24, "man", "z", 2.0e6)
        phi_kelvin = rad.KelvinRadiaScalarPotential(
            coil, tuple(kelvin_center), radius, tuple(physical_center))
        h_kelvin = rad.KelvinRadiaFieldStrength(
            coil, tuple(kelvin_center), radius, tuple(physical_center))

        delta = point - kelvin_center
        mapped = physical_center + radius * radius / np.dot(delta, delta) * delta
        expected = -float(rad.Fld(coil, "phi", [mapped.tolist()])[0])
        assert math.isclose(float(phi_kelvin(mesh(*point))), expected,
                            rel_tol=2.0e-11, abs_tol=2.0e-14)

        # The native coefficient deliberately has no symbolic derivative;
        # assess its differential contract by centered physical evaluation.
        step = 2.0e-5
        numerical_gradient = np.empty(3)
        for component in range(3):
            plus = point.copy()
            minus = point.copy()
            plus[component] += step
            minus[component] -= step
            numerical_gradient[component] = (
                float(phi_kelvin(mesh(*plus))) - float(phi_kelvin(mesh(*minus)))
            ) / (2.0 * step)
        expected_h = np.asarray(h_kelvin(mesh(*point)), dtype=float)
        assert np.linalg.norm(-numerical_gradient - expected_h) / np.linalg.norm(expected_h) < 5.0e-4
    finally:
        rad.UtiDelAll()
