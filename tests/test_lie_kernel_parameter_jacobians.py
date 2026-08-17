"""Lock the no-jacobian mode of the fourth-order Lie kernel.

``parameter_jacobians=False`` must reproduce the map tensors bit-for-bit
(the value path is identical; only the topology-optimization adjoint over
``40 * n_segment`` parameters is skipped) while returning zero-width
jacobian arrays.  The adjoint composition dominates the kernel runtime and
scales superlinearly in the segment count, so verification callers rely on
this mode.
"""

import numpy as np

from radia.accelerator_lie_topopt import (
    _fourth_order_lie_map_from_vector_potential_polynomials,
)

RIGIDITY = 3.0
H = 0.125


def synthetic_arrays(n_seg, degree=5):
    Ay = np.zeros((n_seg, degree + 1, degree + 1))
    As = np.zeros((n_seg, degree + 1, degree + 1))
    As[:, 1, 0] = -H * RIGIDITY          # dipole balancing the curvature
    As[:, 2, 0] = 0.5 * 8.0              # quadrupole-like term
    As[:, 0, 2] = -0.5 * 8.0             # its harmonic y^2 partner
    As[:, 3, 0] = 20.0                   # mild sextupole for f3 content
    As[:, 1, 2] = -3 * 20.0
    lengths = np.full(n_seg, 0.02 / n_seg)
    curvatures = np.full(n_seg, H)
    return Ay, As, lengths, curvatures


def test_no_jacobian_mode_matches_and_zero_width():
    Ay, As, lengths, curvatures = synthetic_arrays(4)
    kwargs = dict(
        reference_curvature_per_m=curvatures,
        longitudinal_component="covariant",
        reference_orbit_tolerance=1.0e-6,
    )
    with_jac = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay, As, lengths, RIGIDITY, **kwargs)
    without = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay, As, lengths, RIGIDITY, parameter_jacobians=False, **kwargs)

    for name in ("R", "T", "U", "V", "f3", "f4", "f5"):
        a = getattr(with_jac.transfer, name)
        b = getattr(without.transfer, name)
        assert np.array_equal(a, b), f"tensor {name} changed"
    assert np.array_equal(with_jac.hamiltonian_linear,
                          without.hamiltonian_linear)
    assert without.transfer.R_jacobian.shape[0] == 0
    assert without.transfer.f3_jacobian.shape[0] == 0
    assert without.hamiltonian_linear_jacobian.shape[0] == 0
    assert with_jac.transfer.R_jacobian.shape[0] == len(
        with_jac.parameter_names) * lengths.size


def test_spoly_constant_matches_constant_path():
    # k_s = 0 s-polynomial arrays through the nonautonomous path must
    # reproduce the constant path (identical stage arithmetic).
    Ay, As, lengths, curvatures = synthetic_arrays(3)
    kwargs = dict(
        reference_curvature_per_m=curvatures,
        longitudinal_component="covariant",
        reference_orbit_tolerance=1.0e-6,
        parameter_jacobians=False,
    )
    constant = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay, As, lengths, RIGIDITY, **kwargs)
    spoly = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay[:, None], As[:, None], lengths, RIGIDITY, **kwargs)
    for name in ("R", "T", "U", "V", "f3", "f4", "f5"):
        a = getattr(constant.transfer, name)
        b = getattr(spoly.transfer, name)
        assert np.allclose(a, b, rtol=0.0, atol=1.0e-14), name


def test_spoly_fourth_order_in_segment_s_dependence():
    # One segment with linear+quadratic zeta-dependence of the quadrupole
    # strength: the nonautonomous flow must converge at fourth order in
    # the step, and beat single-jet midpoint staging on the same segment.
    degree = 5
    length = 0.05

    def arrays(ks):
        Ay = np.zeros((1, ks + 1, degree + 1, degree + 1))
        As = np.zeros((1, ks + 1, degree + 1, degree + 1))
        As[0, 0, 1, 0] = -H * RIGIDITY
        for k, strength in ((0, 8.0), (1, 5.0), (2, -4.0)):
            if k <= ks:
                As[0, k, 2, 0] = 0.5 * strength
                As[0, k, 0, 2] = -0.5 * strength
        return Ay, As

    Ay, As = arrays(2)
    kwargs = dict(
        reference_curvature_per_m=np.array([H]),
        longitudinal_component="covariant",
        reference_orbit_tolerance=1.0e-6,
        parameter_jacobians=False,
        # Deliberately coarse maps measure the flow truncation itself, so
        # the Dragt-Finn consistency gate must not reject them here.
        factorization_tolerance=1.0e-3,
    )

    def spoly_map(step):
        return _fourth_order_lie_map_from_vector_potential_polynomials(
            Ay, As, np.array([length]), RIGIDITY,
            maximum_step_m=step, **kwargs).transfer

    reference = spoly_map(length / 64.0)
    coarse = spoly_map(length / 2.0)
    fine = spoly_map(length / 4.0)
    err_coarse = float(np.max(np.abs(coarse.R - reference.R)))
    err_fine = float(np.max(np.abs(fine.R - reference.R)))
    assert err_fine < err_coarse / 10.0     # ~16x for a clean fourth order

    # Midpoint staging (single constant jet at zeta=0, where the k>=1
    # terms vanish) misses the s-dependence entirely: its error must
    # dwarf the converged spoly result.
    midpoint = _fourth_order_lie_map_from_vector_potential_polynomials(
        Ay[:, 0], As[:, 0], np.array([length]), RIGIDITY,
        maximum_step_m=length / 64.0, **kwargs).transfer
    err_midpoint = float(np.max(np.abs(midpoint.R - reference.R)))
    assert err_midpoint > 30.0 * max(err_fine, 1.0e-300)
