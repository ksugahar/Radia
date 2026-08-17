"""Construction locks for the CanonicalHCurl space (radia.beam_canonical_hcurl).

Locks the measured/oracle-verified structure: the dimension law
``dim = p_x*(p_s+1)``, curved-chart fixed-dimension definition, machine-
precision representation of exact vacuum fields, the graded L1 interface
contract's spline counts, hard interface continuity after a chain fit, and
the midplane-only rejection.
"""

import numpy as np
import pytest

from radia.beam_canonical_hcurl import (
    CanonicalHCurlChain,
    CanonicalHCurlElement,
    graded_breaks,
)

HW, HH = 0.010, 0.0035


def make_element(order_x, order_s, curvature=(0.0,), half_length=0.005):
    return CanonicalHCurlElement(
        order_x=order_x, order_s=order_s,
        half_width_m=HW, half_height_m=HH, half_length_m=half_length,
        curvature_poly_per_m=curvature,
    )


def test_dimension_law_flat_and_curved():
    for order_x, order_s in ((4, 0), (6, 1), (6, 3), (8, 2)):
        for curvature in ((0.0,), (0.125,), (5.0,), (0.1, 0.05)):
            if len(curvature) > order_s + 1:
                continue
            element = make_element(order_x, order_s, curvature)
            assert element.dimension == order_x * (order_s + 1)
            assert element.basis.shape == (
                len(element.ay_exponents) + len(element.as_exponents),
                element.dimension,
            )


def test_curved_defects_are_small_perturbations():
    element = make_element(8, 2, curvature=(0.125,))
    # Retained beyond-strict-kernel defects are O(h * truncation tail)
    # relative to the constraint operator norm (measured ~4e-5 at h=0.125).
    relative = float(np.max(element.vacuum_defects)) / element.vacuum_defect_scale
    assert relative < 1.0e-4
    flat = make_element(8, 2)
    flat_relative = float(np.max(flat.vacuum_defects)) / flat.vacuum_defect_scale
    assert flat_relative < 1.0e-12


def test_uniform_vertical_field_exact_in_curved_chart():
    h = 0.125
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.01]), HW, HH, order_x=4, order_s=1,
        curvature_per_m=lambda s: h)
    rng = np.random.default_rng(7)
    n = 4 * chain.chain_dimension + 200
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(0.0, 0.01, n)
    fit = chain.fit_frame_samples(x, y, s, np.zeros(n),
                                  np.full(n, 0.375), np.zeros(n))
    assert fit.maximum_residual_t < 1.0e-12
    b = chain.magnetic_flux_density_frame(x, y, s)
    assert float(np.max(np.abs(b - [0.0, 0.375, 0.0]))) < 1.0e-12


def test_two_dimensional_multipoles_machine_precision():
    element_chain = CanonicalHCurlChain(
        np.array([-0.005, 0.005]), HW, HH, order_x=6, order_s=0)
    rng = np.random.default_rng(11)
    n = 4 * element_chain.chain_dimension + 400
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = np.zeros(n)
    z = x + 1j * y
    for m in (1, 2, 3, 4):
        # a_s = -Re[z^(m+1)]/(m+1)  =>  (B_x, B_y) = (Im[z^m], Re[z^m]).
        by = np.real(z**m) / HW**m            # normalized amplitude
        bx = np.imag(z**m) / HW**m
        fit = element_chain.fit_frame_samples(x, y, s, bx, by, np.zeros(n))
        assert fit.maximum_residual_t < 1.0e-12, f"multipole m={m}"


def test_chain_dimension_matches_spline_law():
    # Verified through p_s=5: the L1-contract strict kernel at h=0 equals
    # the fixed-dimension target exactly (machine-zero interface defects),
    # so the spline law p_x*(E+p_s) is the measured truth, not a guess.
    for order_s in (2, 3, 4, 5):
        chain = CanonicalHCurlChain(
            np.linspace(0.0, 0.08, 9), HW, HH, order_x=8, order_s=order_s)
        assert chain.element_count == 8
        # Per multipole: E*(p_s+1) - (E-1)*conditions = E + p_s.
        assert chain.chain_dimension == 8 * (8 + order_s)
        if chain.interface_defects.size:
            relative = float(np.max(chain.interface_defects)) \
                / chain.interface_defect_scale
            assert relative < 1.0e-12


def test_chain_dimension_law_survives_varying_curvature():
    # Per-element curvature polynomials break the flat trace alignment, so
    # the strict joint kernel collapses; the fixed-dimension least-defect
    # reduction keeps the spline law with O(h) interface defects.
    chain = CanonicalHCurlChain(
        np.linspace(0.0, 0.08, 9), HW, HH, order_x=8, order_s=2,
        curvature_per_m=lambda s: 0.05 + 1.0 * s)
    assert chain.chain_dimension == 8 * (8 + 2)
    assert chain.interface_defects.size == chain.chain_dimension
    relative = float(np.max(chain.interface_defects)) \
        / chain.interface_defect_scale
    assert relative < 1.0e-2


def test_order_s_one_chain_rejected():
    with pytest.raises(ValueError, match="order_s >= 2"):
        CanonicalHCurlChain(np.linspace(0.0, 0.03, 4), HW, HH,
                            order_x=6, order_s=1)


def test_chain_fit_interfaces_are_hard():
    y0 = 0.02
    scale = 2.0e-7

    def field(x, y, s):
        r2p = x**2 + (y - y0)**2
        r2m = x**2 + (y + y0)**2
        bx = -scale * ((y - y0) / r2p + (y + y0) / r2m)
        by = scale * (x / r2p + x / r2m)
        return bx, by, np.zeros_like(x)

    chain = CanonicalHCurlChain(
        np.linspace(-0.02, 0.02, 5), HW, HH, order_x=8, order_s=2)
    rng = np.random.default_rng(23)
    n = 4 * chain.chain_dimension + 800
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(-0.02, 0.02, n)
    bx, by, bs = field(x, y, s)
    fit = chain.fit_frame_samples(x, y, s, bx, by, bs)
    assert fit.maximum_interface_ay_jump < 1.0e-14
    assert fit.maximum_interface_b_value_jump < 1.0e-14
    assert fit.relative_residual < 1.0e-3
    b = chain.magnetic_flux_density_frame(x, y, s)
    check = np.max(np.abs(b - np.column_stack((bx, by, bs))))
    assert float(check) < 4.0 * fit.maximum_residual_t + 1e-15


def test_sample_starvation_rejected():
    # One thin element with (almost) no samples must fail loud, not fit.
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.001, 0.04]), HW, HH, order_x=6, order_s=2)
    rng = np.random.default_rng(9)
    n = 4 * chain.chain_dimension + 400
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(0.002, 0.04, n)   # never inside the thin first element
    with pytest.raises(ValueError, match="sample starvation"):
        chain.fit_frame_samples(x, y, s, np.zeros(n), np.ones(n),
                                np.zeros(n))


def test_midplane_only_cloud_rejected():
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.01]), HW, HH, order_x=4, order_s=0)
    rng = np.random.default_rng(5)
    n = 4 * chain.chain_dimension + 100
    x = rng.uniform(-HW, HW, n)
    y = np.zeros(n)
    s = rng.uniform(0.0, 0.01, n)
    with pytest.raises(ValueError, match="midplane-only"):
        chain.fit_frame_samples(x, y, s, np.zeros(n), np.ones(n), np.zeros(n))


def test_graded_breaks_equidistribute_monitor():
    s = np.linspace(0.0, 1.0, 401)
    # Two fringe bumps near 0.2 and 0.8.
    w = np.exp(-((s - 0.2) / 0.03) ** 2) + np.exp(-((s - 0.8) / 0.03) ** 2)
    uniform = graded_breaks(s, np.zeros_like(s), 10)
    assert np.allclose(uniform, np.linspace(0.0, 1.0, 11))
    graded = graded_breaks(s, w, 10, strength=8.0)
    assert graded[0] == 0.0 and graded[-1] == 1.0
    sizes = np.diff(graded)
    def local_size(position):
        index = int(np.clip(np.searchsorted(graded, position) - 1, 0, 9))
        return sizes[index]
    # Elements shrink at the bumps relative to the flat middle.
    assert local_size(0.2) < 0.4 * local_size(0.5)
    assert local_size(0.8) < 0.4 * local_size(0.5)


def test_gradient_matches_finite_difference():
    chain = CanonicalHCurlChain(
        np.array([-0.01, 0.01]), HW, HH, order_x=6, order_s=2,
        curvature_per_m=lambda s: 0.1 + 2.0 * s)
    rng = np.random.default_rng(17)
    n = 4 * chain.chain_dimension + 300
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(-0.01, 0.01, n)
    z = x + 1j * y
    chain.fit_frame_samples(x, y, s, np.imag(z**2) / HW**2,
                            np.real(z**2) / HW**2, np.zeros(n))
    px_, py_, ps_ = 0.004, 0.001, 0.003
    a0, grad = chain.vector_potential_and_gradient_frame(
        np.array([px_]), np.array([py_]), np.array([ps_]))
    step = 1.0e-6
    for axis, offset in ((0, (step, 0.0)), (1, (0.0, step))):
        plus = chain.vector_potential_frame(
            np.array([px_ + offset[0]]), np.array([py_ + offset[1]]),
            np.array([ps_]))
        minus = chain.vector_potential_frame(
            np.array([px_ - offset[0]]), np.array([py_ - offset[1]]),
            np.array([ps_]))
        fd = (plus - minus)[0] / (2.0 * step)
        assert float(np.max(np.abs(grad[0, :, axis] - fd))) < 1.0e-8


def test_curl_of_a_matches_b_evaluator_in_curved_chart():
    # Internal consistency lock: the B evaluator IS the exact polynomial
    # curl of the fitted A in the curved chart.  Reconstruct curl A by
    # chart finite differences (Bx=(dAs/dy-dAy/ds)/g, By=-dAs/dx/g,
    # Bs=dAy/dx with covariant As) and compare.
    h = 0.125
    chain = CanonicalHCurlChain(
        np.array([-0.01, 0.01]), HW, HH, order_x=6, order_s=2,
        curvature_per_m=lambda s: h)
    rng = np.random.default_rng(41)
    n = 4 * chain.chain_dimension + 300
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(-0.01, 0.01, n)
    z = x + 1j * y
    chain.fit_frame_samples(x, y, s, np.imag(z**3) / HW**3,
                            np.real(z**3) / HW**3, np.zeros(n))
    probes = np.column_stack((
        rng.uniform(-0.5 * HW, 0.5 * HW, 20),
        rng.uniform(-0.5 * HH, 0.5 * HH, 20),
        rng.uniform(-0.008, 0.008, 20),
    ))
    step = 1.0e-6

    def a_at(px_, py_, ps_):
        return chain.vector_potential_frame(px_, py_, ps_)

    b_eval = chain.magnetic_flux_density_frame(
        probes[:, 0], probes[:, 1], probes[:, 2])
    for row, (px_, py_, ps_) in enumerate(probes):
        d_dx = (a_at([px_ + step], [py_], [ps_])
                - a_at([px_ - step], [py_], [ps_]))[0] / (2 * step)
        d_dy = (a_at([px_], [py_ + step], [ps_])
                - a_at([px_], [py_ - step], [ps_]))[0] / (2 * step)
        d_ds = (a_at([px_], [py_], [ps_ + step])
                - a_at([px_], [py_], [ps_ - step]))[0] / (2 * step)
        g = 1.0 + h * px_
        curl_fd = np.array([
            (d_dy[2] - d_ds[1]) / g,
            -d_dx[2] / g,
            d_dx[1],
        ])
        assert float(np.max(np.abs(curl_fd - b_eval[row]))) < 1.0e-7


def test_periodic_chain_spline_dimension_and_cohomology():
    for order_s in (2, 3, 4, 5):
        ring = CanonicalHCurlChain(
            np.linspace(0.0, 0.08, 9), HW, HH, order_x=6, order_s=order_s,
            periodic=True)
        # Maximal-smoothness periodic spline: E DOFs per multipole family.
        assert ring.chain_dimension == 6 * 8

    ring = CanonicalHCurlChain(
        np.linspace(0.0, 0.08, 9), HW, HH, order_x=6, order_s=2,
        periodic=True)
    rng = np.random.default_rng(31)
    n = 4 * ring.chain_dimension + 400
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(0.0, 0.08, n)
    fit = ring.fit_frame_samples(x, y, s, np.zeros(n), np.full(n, 0.2),
                                 np.zeros(n))
    b_before = ring.magnetic_flux_density_frame(x, y, s)
    a_before = ring.vector_potential_frame(x, y, s)
    ring.set_ring_circulation(1.6e-4)      # linked flux, T*m^2
    b_after = ring.magnetic_flux_density_frame(x, y, s)
    a_after = ring.vector_potential_frame(x, y, s)
    # Cohomology mode: pure circulation, zero field.
    assert float(np.max(np.abs(b_after - b_before))) == 0.0
    assert np.allclose(a_after[:, 2] - a_before[:, 2], 1.6e-4 / 0.08)
    assert float(np.max(np.abs(a_after[:, 1] - a_before[:, 1]))) == 0.0
    assert fit.maximum_residual_t < 1.0e-12

    with pytest.raises(ValueError, match="periodic"):
        open_chain = CanonicalHCurlChain(
            np.linspace(0.0, 0.08, 9), HW, HH, order_x=6, order_s=2)
        open_chain.set_ring_circulation(1.0e-4)


def test_lie_element_spoly_arrays_match_pointwise_coefficients():
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.01, 0.02]), HW, HH, order_x=5, order_s=2,
        curvature_per_m=lambda s: 0.1 + 3.0 * s)
    rng = np.random.default_rng(53)
    n = 4 * chain.chain_dimension + 400
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(0.0, 0.02, n)
    z = x + 1j * y
    chain.fit_frame_samples(x, y, s, np.imag(z**2) / HW**2,
                            np.real(z**2) / HW**2, np.zeros(n))
    Ay, As, lengths, curvature_polys = chain.lie_element_spoly_arrays(5)
    assert Ay.shape == (2, 3, 6, 6) and np.allclose(lengths, 0.01)
    # Collapsing the zeta-polynomial must reproduce the midpoint-staging
    # coefficient extraction at any zeta.
    offsets = [0, chain.elements[0].dimension]
    for e, zeta in ((0, -0.3), (1, 0.7)):
        powers = zeta ** np.arange(3)
        Ay_z = np.tensordot(powers, Ay[e], 1)
        As_z = np.tensordot(powers, As[e], 1)
        Ay_map, As_map = chain.elements[e].transverse_coefficients(zeta, 5)
        c = chain._fit.coefficients[
            offsets[e]:offsets[e] + chain.elements[e].dimension]
        assert np.allclose(Ay_z, Ay_map @ c, atol=1e-15)
        assert np.allclose(As_z, As_map @ c, atol=1e-15)
        h_z = float(powers @ curvature_polys[e])
        expected = 0.1 + 3.0 * (0.005 + 0.01 * e + 0.005 * zeta)
        assert abs(h_z - expected) < 1e-9


def test_lie_segment_arrays_contract():
    h = 0.1
    chain = CanonicalHCurlChain(
        np.array([0.0, 0.02]), HW, HH, order_x=5, order_s=2,
        curvature_per_m=lambda s: h)
    rng = np.random.default_rng(3)
    n = 4 * chain.chain_dimension + 300
    x = rng.uniform(-HW, HW, n)
    y = rng.uniform(-HH, HH, n)
    s = rng.uniform(0.0, 0.02, n)
    chain.fit_frame_samples(x, y, s, np.zeros(n), np.full(n, 0.2),
                            np.zeros(n))
    Ay, As, lengths, curvatures = chain.lie_segment_arrays(4, degree=5)
    assert Ay.shape == (4, 6, 6) and As.shape == (4, 6, 6)
    assert np.allclose(lengths, 0.005)
    assert np.allclose(curvatures, h)
    # Design-orbit gauge zeros are structural.
    assert float(np.max(np.abs(Ay[:, 0, 0]))) == 0.0
    assert float(np.max(np.abs(As[:, 0, 0]))) == 0.0
    # Coefficient arrays reproduce the evaluator (flat metric direction x=0
    # so covariant == physical for the As comparison at x=0).
    mids = np.array([0.0025 + 0.005 * k for k in range(4)])
    probe_y = 0.001
    a_eval = chain.vector_potential_frame(
        np.zeros(4), np.full(4, probe_y), mids)
    powers = probe_y ** np.arange(6)
    ay_from_coef = Ay[:, 0, :] @ powers
    as_from_coef = As[:, 0, :] @ powers
    assert float(np.max(np.abs(a_eval[:, 1] - ay_from_coef))) < 1.0e-12
    assert float(np.max(np.abs(a_eval[:, 2] - as_from_coef))) < 1.0e-12
