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
