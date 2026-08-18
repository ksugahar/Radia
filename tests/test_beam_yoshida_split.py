"""Locks for the Yoshida split-operator tracker on CanonicalHCurl chains.

The integrator's contract: exact symplectic factors at ANY amplitude and
step (machine-level M^T J M = J), the deck's dipole balance (on-orbit
state stays on orbit), order-2/4 convergence in the step count, and
agreement with the exact-sqrt canonical A-RK at small amplitude where the
paraxial gap is negligible.
"""

import numpy as np
import pytest
from scipy.integrate import solve_ivp

from radia.accelerator_lie_topopt import (
    canonical_vector_potential_hamiltonian_rhs,
)
from radia.beam_canonical_hcurl import CanonicalHCurlChain, CanonicalHCurlFit
from radia.beam_yoshida_split import track_yoshida_split

HW, HH = 0.010, 0.004
RIGIDITY = 3.0


def synthetic_chain(seed=3, curvature=lambda s: 0.12 + 1.5 * s,
                    amplitude=1.0e-3):
    """A legitimate nontrivial chain: random coefficients in the REDUCED
    space (interface contract satisfied, each element exactly in its
    vacuum space), populated s-dependence so ``a_y != 0``."""
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.02, 0.04]), HW, HH, order_x=5, order_s=2,
        curvature_per_m=curvature)
    rng = np.random.default_rng(seed)
    reduced = amplitude * rng.standard_normal(chain.chain_dimension)
    coefficients = chain._reduced @ reduced
    chain._fit = CanonicalHCurlFit(
        coefficients=coefficients, maximum_residual_t=0.0, field_scale_t=1.0,
        maximum_interface_ay_jump=0.0, maximum_interface_b_value_jump=0.0,
        sample_count=0)
    return chain


def test_symplectic_at_large_amplitude_and_single_big_step():
    chain = synthetic_chain()
    state = np.array([5.0e-3, 5.0e-2, 2.0e-3, 3.0e-2])   # LARGE amplitude
    span = (0.0, 0.04)

    def push(z):
        return track_yoshida_split(chain, RIGIDITY, z, span, step_count=2)

    step = 1.0e-6
    jacobian = np.empty((4, 4))
    for column in range(4):
        plus = state.copy()
        minus = state.copy()
        plus[column] += step
        minus[column] -= step
        jacobian[:, column] = (push(plus) - push(minus)) / (2.0 * step)
    J = np.array([[0.0, 1.0, 0.0, 0.0], [-1.0, 0.0, 0.0, 0.0],
                  [0.0, 0.0, 0.0, 1.0], [0.0, 0.0, -1.0, 0.0]])
    residual = jacobian.T @ J @ jacobian - J
    assert float(np.max(np.abs(residual))) < 5.0e-7   # FD-limited


def test_dipole_balance_keeps_the_orbit():
    by = 0.375
    htilde = by / RIGIDITY
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.02, 0.04]), HW, HH, order_x=4, order_s=2,
        curvature_per_m=lambda s: htilde)
    rng = np.random.default_rng(7)
    n = 4 * chain.chain_dimension + 300
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(0.0, HH, n)
    s = rng.uniform(0.0, 0.04, n)
    chain.fit_frame_samples(x, y, s, np.zeros(n), np.full(n, by), np.zeros(n))
    out = track_yoshida_split(chain, RIGIDITY, np.zeros(4), (0.0, 0.04),
                              step_count=16)
    assert float(np.max(np.abs(out))) < 1.0e-12


def test_step_convergence_orders_two_and_four():
    chain = synthetic_chain()
    state = np.array([2.0e-3, 2.0e-3, 1.0e-3, 1.5e-3])
    span = (0.0, 0.04)
    reference = track_yoshida_split(chain, RIGIDITY, state, span,
                                    step_count=256, order=4)

    def error(steps, order):
        value = track_yoshida_split(chain, RIGIDITY, state, span,
                                    step_count=steps, order=order)
        return float(np.max(np.abs(value - reference)))

    second_coarse = error(8, 2)
    second_fine = error(16, 2)
    assert second_coarse / second_fine > 3.0      # ~4x for clean order 2
    fourth_coarse = error(4, 4)
    fourth_fine = error(8, 4)
    assert fourth_coarse / fourth_fine > 10.0     # ~16x for clean order 4


def balanced_dipole_chain():
    """Uniform By with the matching curvature: the reference orbit is a
    true trajectory, so the paraxial gap is driven by the INITIAL
    momentum alone (a synthetic chain without the dipole balance pumps
    p to ~integral(htilde) regardless of field and buries the scaling)."""
    by = 0.375
    htilde = by / RIGIDITY
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.02, 0.04]), HW, HH, order_x=4, order_s=2,
        curvature_per_m=lambda s: htilde)
    rng = np.random.default_rng(7)
    n = 4 * chain.chain_dimension + 300
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(0.0, HH, n)
    s = rng.uniform(0.0, 0.04, n)
    chain.fit_frame_samples(x, y, s, np.zeros(n), np.full(n, by), np.zeros(n))
    return chain


def paraxial_gap(chain, momentum):
    """Converged order-4 split vs the exact-sqrt canonical A-RK."""
    span = (0.0, 0.04)
    state = np.array([0.0, momentum, 0.0, 0.4 * momentum])

    def rhs(s_value, z):
        a, gradient = chain.vector_potential_and_gradient_frame(
            np.array([float(z[0])]), np.array([float(z[2])]),
            np.array([s_value]))
        htilde = float(chain.curvature_frame(np.array([s_value]))[0])
        return canonical_vector_potential_hamiltonian_rhs(
            z, a[0] / RIGIDITY, gradient[0] / RIGIDITY,
            reference_curvature_per_m=htilde)

    six = np.array([state[0], state[1], state[2], state[3], 0.0, 0.0])
    solution = solve_ivp(rhs, span, six, method="DOP853",
                         rtol=1.0e-12, atol=1.0e-13)
    assert solution.success
    a_rk = solution.y[:4, -1]
    split = track_yoshida_split(chain, RIGIDITY, state, span,
                                step_count=128, order=4)
    return float(np.max(np.abs(split - a_rk)))


def test_paraxial_gap_is_zero_on_orbit_and_scales_as_momentum_squared():
    # The first clean quantification of the deck's paraxial approximation
    # (measured 2026-08-18): zero on the balanced orbit at zero momentum
    # (every dropped term carries p^2), and EXACTLY quadratic in the
    # launch momentum -- ratios 4.000/4.000/4.001 across three doublings,
    # gap ~ 2.9e-7 at 10 mrad over this 80 mm dipole.
    chain = balanced_dipole_chain()
    assert paraxial_gap(chain, 0.0) < 1.0e-15
    gap_5 = paraxial_gap(chain, 5.0e-3)
    gap_10 = paraxial_gap(chain, 1.0e-2)
    assert 3.6 < gap_10 / gap_5 < 4.4
    assert 2.0e-7 < gap_10 < 4.0e-7


def test_input_validation():
    chain = synthetic_chain()
    with pytest.raises(ValueError, match="order"):
        track_yoshida_split(chain, RIGIDITY, np.zeros(4), (0.0, 0.01),
                            step_count=4, order=3)
    with pytest.raises(ValueError, match="advance"):
        track_yoshida_split(chain, RIGIDITY, np.zeros(4), (0.01, 0.0),
                            step_count=4)
