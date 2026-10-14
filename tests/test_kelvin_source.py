"""Sanity tests for kelvin_source.biot_savart_A_at_points + Kelvin map.

Compares Biot-Savart A_z at the axis of a single circular loop
with the analytical formula (on axis only, since this is where the
closed-form is simple):

    A_theta(rho=R, z=z0) analytical via elliptic integrals -- we
    simplify by testing A at the CENTER of a single small loop which
    is zero by symmetry except for the self-induced component.

More useful sanity check: compare two different discretizations of
the same coil for self-consistency.
"""

import math

import numpy as np

MU_0 = 4e-7 * math.pi


def build_loop(R, z0, n, current):
    """Discretize a circular loop in xy-plane at z=z0 into n segments."""
    thetas = np.linspace(0, 2 * math.pi, n + 1)
    pts = np.column_stack([
        R * np.cos(thetas), R * np.sin(thetas),
        np.full_like(thetas, z0)])
    segs = [(pts[i].tolist(), pts[i + 1].tolist()) for i in range(n)]
    return [segs], [current]


def test_on_axis_B_circular_loop():
    """Analytical B_z on axis of a single circular loop:

        B_z(0, 0, z) = mu_0 * I * R^2 / (2 * (R^2 + z^2)^{3/2})

    We compute A at (eps, 0, z) and take finite-difference curl to
    get B_z ~ (A_theta_at_z+dz - A_theta_at_z-dz) / (2 dz) -- too
    fragile. Instead: just check that A_phi scales ~ R^2/(R^2+z^2)^{1/2}
    magnitude.
    """
    from radia.kelvin_source import biot_savart_A_at_points

    R = 0.010
    I = 1.0
    paths, currents = build_loop(R, 0.0, 200, I)
    # Evaluate A on axis at several z
    for z in [0.0, 0.005, 0.010, 0.020]:
        obs = np.array([[1e-6, 0.0, z]])
        A = biot_savart_A_at_points(obs, paths, currents, n_quad=6)
        # On axis, dipole expansion for small rho:
        #   A_phi(rho, z) ~ mu_0 m * rho / (4 pi (R^2 + z^2)^{3/2})
        # where m = I * pi * R^2 (magnetic moment), so
        #   A_phi ~ (mu_0 I / 4) * rho * R^2 / (R^2 + z^2)^{3/2}
        rho = 1e-6
        expected_A_y = (MU_0 * I / 4.0) * rho * R ** 2 \
            / (R ** 2 + z ** 2) ** 1.5
        A_y = A[0, 1]
        rel_err = abs(A_y - expected_A_y) / abs(expected_A_y)
        print(f"  z={z*1e3:5.1f}mm   A_y={A_y:.4e}   "
              f"expected={expected_A_y:.4e}   rel_err={rel_err*100:.3f}%")
        assert rel_err < 0.02, f"Loop A at z={z} off by {rel_err*100:.2f}%"


def test_kelvin_map_involution():
    from radia.kelvin_source import kelvin_map_3d

    center = np.array([0.0, 0.0, 0.0])
    R = 1.0
    rng = np.random.default_rng(42)
    pts = rng.normal(size=(10, 3))
    pts[:, :] *= 2 + np.linalg.norm(pts, axis=1)[:, None]
    mapped = kelvin_map_3d(pts, center, R)
    double = kelvin_map_3d(mapped, center, R)
    err = np.max(np.abs(double - pts))
    print(f"  max involution error: {err:.2e}")
    assert err < 1e-10


def test_kelvin_map_on_sphere():
    """Points on the Kelvin sphere map to themselves."""
    from radia.kelvin_source import kelvin_map_3d

    center = np.array([0.1, -0.2, 0.3])
    R = 0.5
    thetas = np.linspace(0, 2 * math.pi, 12, endpoint=False)
    pts = center + R * np.column_stack([np.cos(thetas),
                                         np.sin(thetas),
                                         np.zeros_like(thetas)])
    mapped = kelvin_map_3d(pts, center, R)
    err = np.max(np.abs(mapped - pts))
    print(f"  on-sphere invariance error: {err:.2e}")
    assert err < 1e-12


def test_scalar_factor_at_sphere_boundary():
    """Kelvin factor (R/rho')^2 = 1 at rho' = R."""
    from radia.kelvin_source import kelvin_factor_scalar

    center = np.array([0.0, 0.0, 0.0])
    R = 0.5
    pts = np.array([[R, 0, 0], [0, R, 0], [0, 0, R], [-R, 0, 0]])
    f = kelvin_factor_scalar(pts, center, R)
    print(f"  factors at rho'=R: {f}")
    assert np.allclose(f, 1.0)


def test_A_s_kelvin_pathway():
    """A_s_at_obs_with_kelvin: inner points unchanged, kelvin_ext points mapped."""
    from radia.kelvin_source import A_s_at_obs_with_kelvin, biot_savart_A_at_points

    R_loop = 0.030
    R_kelvin = 0.060
    offset = np.array([0.150, 0.0, 0.0])   # Kelvin sphere center (Sugahara 2022)
    I = 1.0
    paths, currents = build_loop(R_loop, 0.0, 200, I)
    # Inner point (physical inner domain, far from Kelvin sphere)
    p_in = np.array([0.010, 0.0, 0.005])
    # Kelvin exterior domain point (computational frame): offset +/- dx within radius
    p_out_kelvin_ext = offset + np.array([0.020, 0.005, 0.010])
    # Its Kelvin-inverse physical position:
    from radia.kelvin_source import kelvin_map_3d
    p_out_phys = kelvin_map_3d(p_out_kelvin_ext, offset, R_kelvin)

    obs = np.vstack([p_in, p_out_kelvin_ext])
    A_full = A_s_at_obs_with_kelvin(
        obs, paths, currents,
        kelvin_center=offset.tolist(), R_kelvin=R_kelvin,
        factor_mode='scalar')

    # Inner vs direct
    A_in_direct = biot_savart_A_at_points(p_in.reshape(1, 3),
                                           paths, currents, n_quad=8)[0]
    err_in = np.max(np.abs(A_full[0] - A_in_direct))
    print(f"  inner A match: err = {err_in:.2e}")
    assert err_in < 1e-12

    # Kelvin exterior domain: expected = A(physical_equivalent) * (R/rho')^2
    A_at_phys = biot_savart_A_at_points(p_out_phys.reshape(1, 3),
                                         paths, currents, n_quad=8)[0]
    rho_prime = np.linalg.norm(p_out_kelvin_ext - offset)
    factor = (R_kelvin / rho_prime) ** 2
    expected = A_at_phys * factor
    err_sh = np.max(np.abs(A_full[1] - expected))
    rel_sh = err_sh / np.max(np.abs(expected))
    print(f"  kelvin_ext factor match: err = {err_sh:.2e}  rel = {rel_sh:.2e}")
    # Quadrature-level tolerance (loop discretization + n_quad=8).
    assert rel_sh < 1e-6


def test_pullback_equals_scalar_for_uniform_A():
    """For A with no radial component, pullback == scalar * A (uniform case)."""
    from radia.kelvin_source import kelvin_factor_scalar, kelvin_pullback_vector

    center = np.array([0.15, 0.0, 0.0])
    R = 0.060
    # Kelvin exterior domain points
    rng = np.random.default_rng(0)
    r_prime = center + 0.5 * R * rng.standard_normal((20, 3))
    # A with no radial component: A perpendicular to (r' - c).
    # For a uniform B = B_z z_hat background, A = (B_z/2) (-y, x, 0) in 3D
    # is azimuthal and is always perpendicular to the radial direction
    # from z-axis but NOT to the radial from a generic center c.
    # Construct A purely tangential to (r' - c):
    d = r_prime - center
    rho = np.linalg.norm(d, axis=1, keepdims=True)
    n = d / rho
    v = rng.standard_normal((20, 3))
    A_tan = v - (np.sum(v * n, axis=1, keepdims=True)) * n   # remove radial
    # Pullback vs scalar * A_tan
    A_pull = kelvin_pullback_vector(A_tan, r_prime, center, R)
    f = kelvin_factor_scalar(r_prime, center, R)
    A_scalar = f[:, None] * A_tan
    err = np.max(np.abs(A_pull - A_scalar)) / np.max(np.abs(A_scalar))
    print(f"  tangential-A: pullback == scalar factor, rel err = {err:.2e}")
    assert err < 1e-12


def test_pullback_flips_radial_sign():
    """Pure radial A: pullback flips sign and applies (R/rho')^2."""
    from radia.kelvin_source import kelvin_factor_scalar, kelvin_pullback_vector

    center = np.array([0.15, 0.0, 0.0])
    R = 0.060
    rng = np.random.default_rng(1)
    r_prime = center + 0.5 * R * rng.standard_normal((15, 3))
    d = r_prime - center
    rho = np.linalg.norm(d, axis=1, keepdims=True)
    n = d / rho
    # Purely radial A
    mag = rng.standard_normal((15, 1))
    A_rad = mag * n
    A_pull = kelvin_pullback_vector(A_rad, r_prime, center, R)
    f = kelvin_factor_scalar(r_prime, center, R)
    # Expected: -f * A_rad (Householder flips sign of radial component)
    expected = -f[:, None] * A_rad
    err = np.max(np.abs(A_pull - expected)) / np.max(np.abs(expected))
    print(f"  radial-A: pullback == -scalar, rel err = {err:.2e}")
    assert err < 1e-12


def test_pullback_involution():
    """Pulling back twice (with same center, R) returns to original vector.

    Since the Kelvin map is involutive, A''(r) = A(r). Verify numerically.
    """
    from radia.kelvin_source import kelvin_map_3d, kelvin_pullback_vector

    center = np.array([0.0, 0.0, 0.0])
    R = 1.0
    rng = np.random.default_rng(2)
    # Physical points outside the sphere (|r| > R)
    r_phys = rng.standard_normal((10, 3))
    r_phys *= 2 + np.linalg.norm(r_phys, axis=1, keepdims=True)
    A_phys = rng.standard_normal((10, 3))
    # First pullback: A_phys at r_phys -> A' at r' = Kelvin_inv(r_phys)
    r_prime = kelvin_map_3d(r_phys, center, R)
    A_kelvin_ext = kelvin_pullback_vector(A_phys, r_prime, center, R)
    # Second pullback: A_kelvin_ext at r_prime (treated as 'physical') -> back at r
    # For involution to hold, we need to interpret the second pullback as
    # going from the kelvin_ext back to physical. The Jacobian is the same
    # form because phi is self-inverse:
    r_back = kelvin_map_3d(r_prime, center, R)
    err_rback = np.max(np.abs(r_back - r_phys))
    assert err_rback < 1e-10
    A_back = kelvin_pullback_vector(A_kelvin_ext, r_back, center, R)
    err = np.max(np.abs(A_back - A_phys)) / np.max(np.abs(A_phys))
    print(f"  double pullback involution, rel err = {err:.2e}")
    assert err < 1e-12


def test_B_pullback_uniform_field():
    """For uniform B_phys = B_0 z_hat, the pseudovector pullback at any
    r' should give B_comp = -(R/rho')^4 * H * B_phys.

    At points on a coordinate axis, H acts simply (radial flip only),
    so we can check both magnitude and direction in closed form.
    """
    from radia.kelvin_source import kelvin_pullback_B_pseudovector

    center = np.array([0.0, 0.0, 0.0])
    R = 1.0
    B0 = 2.5
    B_phys = np.array([0.0, 0.0, B0])  # uniform B_z

    # Probe along x-axis at various rho': n = (1,0,0), so H flips x;
    # B_phys is in z so H leaves it unchanged.
    for rho in [0.1, 0.3, 0.5, 0.9]:
        rp = np.array([rho, 0.0, 0.0])
        B_comp = kelvin_pullback_B_pseudovector(B_phys, rp, center, R)
        expected = -((R / rho) ** 4) * B_phys
        err = np.max(np.abs(B_comp - expected))
        print(f"  rho={rho}: B_comp={B_comp}, expected={expected}, "
              f"err={err:.2e}")
        assert err < 1e-12


def test_B_pullback_radial_flip():
    """For B_phys purely radial along z-axis (B_phys = B_0 z_hat at r'
    on z-axis), Householder reflects: B_phys.n = B_0, refl = -B_phys.

    So B_comp at z-axis with same B_phys = (-1) * (-B_0 z_hat) * (R/rho)^4
                                       = +B_0 (R/rho)^4 z_hat
    (the -1 sign from -det(H) and the radial flip cancel).
    """
    from radia.kelvin_source import kelvin_pullback_B_pseudovector

    center = np.array([0.0, 0.0, 0.0])
    R = 1.0
    B0 = 1.7
    B_phys = np.array([0.0, 0.0, B0])

    for rho in [0.2, 0.4, 0.6, 0.8]:
        rp = np.array([0.0, 0.0, rho])  # on z-axis
        B_comp = kelvin_pullback_B_pseudovector(B_phys, rp, center, R)
        # n = (0,0,1), B_phys.n = B0, refl = B_phys - 2*B0*z_hat = -B_phys
        # B_comp = -(R/rho)^4 * (-B_phys) = +(R/rho)^4 * B_phys
        expected = ((R / rho) ** 4) * B_phys
        err = np.max(np.abs(B_comp - expected))
        print(f"  rho={rho}: B_comp={B_comp}, expected={expected}, "
              f"err={err:.2e}")
        assert err < 1e-12


def test_curl_A_comp_matches_B_pullback():
    """Internal consistency: curl(A_comp) computed numerically in the
    computational frame should equal the 2-form pullback formula for B.

    Use uniform B_phys = B_0 z_hat (so A_phys = (B_0/2)(-y, x, 0)
    azimuthal, tangential everywhere). Compute A_comp(r') analytically
    via the closed form, then take a 4th-order central finite difference
    for curl_z at a probe point. Compare with the 2-form formula
    B_comp = -(R/rho')^4 * H * B_phys (Householder identity here).
    """
    from radia.kelvin_source import (
        kelvin_pullback_B_pseudovector,
        kelvin_pullback_vector,
    )

    center = np.array([0.0, 0.0, 0.0])
    R = 1.0
    B0 = 1.0

    def A_phys(r):
        # uniform B_z -> azimuthal A
        return np.array([-0.5 * B0 * r[1], 0.5 * B0 * r[0], 0.0])

    def kelvin_inv(rp, c=center, R=R):
        d = rp - c
        rho_sq = float(np.dot(d, d))
        return c + (R * R / rho_sq) * d

    def A_comp(rp):
        rphys = kelvin_inv(rp)
        a = A_phys(rphys)
        return kelvin_pullback_vector(a, rp, center, R)

    # Probe on x-axis at rho' = 0.5
    rp = np.array([0.5, 0.0, 0.0])
    h = 1e-5
    # curl_z = dA_y/dx - dA_x/dy
    dAy_dx = (A_comp(rp + np.array([h, 0, 0]))[1]
              - A_comp(rp - np.array([h, 0, 0]))[1]) / (2 * h)
    dAx_dy = (A_comp(rp + np.array([0, h, 0]))[0]
              - A_comp(rp - np.array([0, h, 0]))[0]) / (2 * h)
    curl_z_numeric = dAy_dx - dAx_dy

    B_comp = kelvin_pullback_B_pseudovector(
        np.array([0.0, 0.0, B0]), rp, center, R)
    print(f"  numerical curl_z(A_comp) = {curl_z_numeric:.6e}")
    print(f"  formula  B_comp_z         = {B_comp[2]:.6e}")
    rel = abs(curl_z_numeric - B_comp[2]) / abs(B_comp[2])
    print(f"  rel err = {rel:.2e}")
    assert rel < 1e-5


def _build_sphere_mesh_centered(center, R, maxh=0.25):
    """Tiny Netgen sphere mesh used only to exercise eval_*_physical_from_gf."""
    from netgen.occ import OCCGeometry, Pnt, Sphere
    from ngsolve import Mesh
    ball = Sphere(Pnt(*center), R)
    ball.mat("outer")
    geo = OCCGeometry(ball)
    ngmesh = geo.GenerateMesh(maxh=maxh)
    return Mesh(ngmesh)


def test_eval_Omega_physical_from_gf():
    """0-form helper: the eval just composes with kelvin_map_3d and reads
    the GridFunction. To exercise only that pipeline (and not FEM
    interpolation error), store a LINEAR CF in computational coords r';
    the helper must return that linear CF evaluated at r' = Kelvin(r_phys).
    """
    try:
        from ngsolve import H1, GridFunction, TaskManager, x, y, z
    except Exception as e:
        print(f"  skipping (no NGSolve): {e}")
        return
    from radia.kelvin_source import eval_Omega_physical_from_gf, kelvin_map_3d

    center = np.array([0.3, -0.2, 0.1])
    R = 1.0
    mesh = _build_sphere_mesh_centered(center, R, maxh=0.25)

    # Omega_comp(r') = r'_x + 2 r'_y + 3 r'_z  (linear in r')
    fes = H1(mesh, order=2)
    gf = GridFunction(fes)
    with TaskManager():
        gf.Set(x + 2.0 * y + 3.0 * z)

        rel_errs = []
        for r_phys_off in [(1.5, 0, 0), (0, 2.0, 0.5), (-1.0, 1.2, -0.8)]:
            r_phys = np.array(r_phys_off) + center
            assert np.linalg.norm(r_phys - center) > R
            r_prime = kelvin_map_3d(r_phys, center, R)
            expected = r_prime[0] + 2.0 * r_prime[1] + 3.0 * r_prime[2]
            got = eval_Omega_physical_from_gf(gf, mesh, r_phys, center, R)
            rel = abs(got - expected) / max(1.0, abs(expected))
            rel_errs.append(rel)
            print(f"  r_phys={r_phys}, r_prime={r_prime}, got={got:.6f}, "
                  f"expected={expected:.6f}, rel={rel:.2e}")
        assert max(rel_errs) < 1e-3, f"0-form eval rel err {max(rel_errs):.2e}"


def test_eval_A_physical_from_gf():
    """1-form helper: store a LINEAR vector A_comp(r') in computational
    coords, then eval_A_physical must return the Kelvin pullback of that
    value at r' = Kelvin(r_phys).
    """
    try:
        from ngsolve import CoefficientFunction, GridFunction, VectorH1, x, y, z
    except Exception as e:
        print(f"  skipping (no NGSolve): {e}")
        return
    from radia.kelvin_source import (
        eval_A_physical_from_gf,
        kelvin_map_3d,
        kelvin_pullback_vector,
    )

    center = np.array([0.2, 0.0, -0.1])
    R = 1.0
    mesh = _build_sphere_mesh_centered(center, R, maxh=0.25)

    # A_comp(r') = linear vector field in r'
    A_cf = CoefficientFunction((y - center[1], -(x - center[0]), z - center[2]))
    fes = VectorH1(mesh, order=2)
    gf = GridFunction(fes)
    gf.Set(A_cf)

    rel_errs = []
    for r_phys_off in [(1.8, 0.0, 0.0), (0.0, -2.1, 0.4), (-1.2, 0.7, -1.0)]:
        r_phys = np.array(r_phys_off) + center
        assert np.linalg.norm(r_phys - center) > R
        r_prime = kelvin_map_3d(r_phys, center, R)
        A_comp_at_rprime = np.array([r_prime[1] - center[1],
                                      -(r_prime[0] - center[0]),
                                      r_prime[2] - center[2]])
        # helper should return the 1-form pullback of A_comp_at_rprime at r_phys
        expected = kelvin_pullback_vector(A_comp_at_rprime, r_phys, center, R)
        got = eval_A_physical_from_gf(gf, mesh, r_phys, center, R)
        rel = np.linalg.norm(got - expected) / max(1.0, np.linalg.norm(expected))
        rel_errs.append(rel)
        print(f"  r_phys={r_phys}, got={got}, expected={expected}, rel={rel:.2e}")
    assert max(rel_errs) < 5e-3, f"1-form eval rel err {max(rel_errs):.2e}"


def test_eval_B_physical_from_gf():
    """2-form helper: store a LINEAR B_comp(r') in computational coords,
    then eval_B_physical must return the 2-form pullback of that value
    at r' = Kelvin(r_phys).
    """
    try:
        from ngsolve import (
            CoefficientFunction,
            GridFunction,
            VectorH1,
            x,
            y,
        )
    except Exception as e:
        print(f"  skipping (no NGSolve): {e}")
        return
    from radia.kelvin_source import (
        eval_B_physical_from_gf,
        kelvin_map_3d,
        kelvin_pullback_B_pseudovector,
    )

    center = np.array([0.0, 0.0, 0.5])
    R = 1.0
    mesh = _build_sphere_mesh_centered(center, R, maxh=0.25)

    B_cf = CoefficientFunction((2.0 * (x - center[0]),
                                  y - center[1] + 0.3,
                                  0.7))
    fes = VectorH1(mesh, order=2)
    gfB = GridFunction(fes)
    gfB.Set(B_cf)

    rel_errs = []
    for r_phys_off in [(1.5, 0.0, 0.0), (0.0, 2.2, 0.3), (-0.8, 0.9, -1.4)]:
        r_phys = np.array(r_phys_off) + center
        assert np.linalg.norm(r_phys - center) > R
        r_prime = kelvin_map_3d(r_phys, center, R)
        B_comp_at_rprime = np.array([2.0 * (r_prime[0] - center[0]),
                                      r_prime[1] - center[1] + 0.3,
                                      0.7])
        expected = kelvin_pullback_B_pseudovector(B_comp_at_rprime, r_phys,
                                                    center, R)
        got = eval_B_physical_from_gf(gfB, mesh, r_phys, center, R)
        rel = np.linalg.norm(got - expected) / max(1.0, np.linalg.norm(expected))
        rel_errs.append(rel)
        print(f"  r_phys={r_phys}, got={got}, expected={expected}, rel={rel:.2e}")
    assert max(rel_errs) < 5e-3, f"2-form eval rel err {max(rel_errs):.2e}"


if __name__ == "__main__":
    print("test_on_axis_B_circular_loop"); test_on_axis_B_circular_loop()
    print("test_kelvin_map_involution"); test_kelvin_map_involution()
    print("test_kelvin_map_on_sphere"); test_kelvin_map_on_sphere()
    print("test_scalar_factor_at_sphere_boundary");
    test_scalar_factor_at_sphere_boundary()
    print("test_A_s_kelvin_pathway"); test_A_s_kelvin_pathway()
    print("test_pullback_equals_scalar_for_uniform_A")
    test_pullback_equals_scalar_for_uniform_A()
    print("test_pullback_flips_radial_sign")
    test_pullback_flips_radial_sign()
    print("test_pullback_involution")
    test_pullback_involution()
    print("test_B_pullback_uniform_field"); test_B_pullback_uniform_field()
    print("test_B_pullback_radial_flip"); test_B_pullback_radial_flip()
    print("test_curl_A_comp_matches_B_pullback")
    test_curl_A_comp_matches_B_pullback()
    print("test_eval_Omega_physical_from_gf")
    test_eval_Omega_physical_from_gf()
    print("test_eval_A_physical_from_gf")
    test_eval_A_physical_from_gf()
    print("test_eval_B_physical_from_gf")
    test_eval_B_physical_from_gf()
    print("\nall tests passed.")
