"""Contract for the matched patch that repairs a surface impedance.

The patch is only worth having if the interior solve behind it is right and if
it refuses to be matched where the outer model is already failing.  Both are
checked here against things that cannot drift: the exact Bessel round wire,
and the geometry of the fixture section.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

MU0 = 4e-7 * math.pi


def _round_wire_panels(a, n=256):
    theta = np.linspace(0.0, 2.0 * np.pi, n, endpoint=False)
    return (np.stack([a * np.cos(theta), a * np.sin(theta)], axis=1),
            np.full(n, 2.0 * np.pi * a / n),
            np.full(n, 1.0 / n, dtype=complex))


def test_interior_solve_reproduces_the_exact_round_wire():
    """The one anchor that is not another solver: the Bessel series.

    The exterior of a round wire is that of a line current at any skin depth,
    so the open boundary is known in closed form; the solve is affine in the
    uniform axial field, so two runs fix the one carrying exactly 1 A.  What
    remains between the solver and the exact resistance is the solver.
    """
    ng = pytest.importorskip("ngsolve")
    jv = pytest.importorskip("scipy.special").jv
    from netgen.occ import WorkPlane

    from radia.sibc_corner_patch import solve_cut_patch

    a, sigma, frequency = 1.0e-3, 5.8e7, 200_000.0
    omega = 2.0 * math.pi * frequency
    delta = math.sqrt(2.0 / (omega * MU0 * sigma))
    panel_xy, panel_ds, panel_current = _round_wire_panels(a)
    radius = ng.sqrt(ng.x * ng.x + ng.y * ng.y)
    common = dict(
        x_cut=None, panel_xy=panel_xy, panel_ds=panel_ds,
        panel_current=panel_current,
        surface_impedance=(1 + 1j) / (sigma * delta), omega=omega,
        sigma=sigma, maxh_conductor=0.4 * delta, air_pad=6.0 * a,
        maxh_air=a, order=2,
        exterior_potential=(MU0 / (2.0 * math.pi)) * ng.log(a / radius))
    face = WorkPlane().Circle(0.0, 0.0, a).Face()

    with pytest.warns(RuntimeWarning, match="unverified boundary data"):
        currents = [solve_cut_patch(face, axial_E_field=drive,
                                    **common).total_current()
                    for drive in (0.0 + 0j, 1.0 + 0j)]
        drive = (1.0 - currents[0]) / (currents[1] - currents[0])
        patch = solve_cut_patch(face, axial_E_field=drive, **common)
    assert patch.boundary_data_error is None

    assert abs(patch.total_current()) == pytest.approx(1.0, rel=1e-9)
    resistance = 2.0 * patch.total_loss() / abs(patch.total_current()) ** 2
    k = np.sqrt(-1j * omega * MU0 * sigma)
    exact = float(((k / (2.0 * math.pi * a * sigma))
                   * jv(0, k * a) / jv(1, k * a)).real)
    assert resistance == pytest.approx(exact, rel=2e-3)


def test_the_patch_refuses_a_cut_through_the_fin():
    """A match made inside the taper is a match where the outer model fails.

    The whole construction rests on the cut chord carrying the flat half-space
    skin profile, which is a statement about the thick body.  Cutting through
    the beak would take the outer model's word exactly where it is known to be
    wrong, so it is refused rather than silently allowed.
    """
    pytest.importorskip("ngsolve")
    from netgen.occ import MoveTo

    from radia.sibc_corner_patch import solve_cut_patch

    # A wedge: full height at x = 0, tapering to a point at x = 4 mm.
    outline = np.array([[0.0, -2.0], [4.0, 0.0], [0.0, 2.0]]) * 1e-3
    dense = np.concatenate([
        np.linspace(outline[i], outline[(i + 1) % 3], 64, endpoint=False)
        for i in range(3)])
    ds = np.linalg.norm(dense - np.roll(dense, 1, axis=0), axis=1)
    face = MoveTo(0.0, -2e-3).Rectangle(4e-3, 4e-3).Face()

    with pytest.raises(ValueError, match="full-thickness body"):
        solve_cut_patch(
            face, x_cut=2.0e-3, panel_xy=dense, panel_ds=ds,
            panel_current=np.full(len(dense), 1.0 / len(dense), dtype=complex),
            axial_E_field=1.0 + 0j, surface_impedance=1e-3 + 1e-3j,
            omega=2.0 * math.pi * 150_000.0, sigma=5.8e7,
            maxh_conductor=1e-4)


def test_the_patch_refuses_a_body_thinner_than_its_own_skin():
    """Matching needs the outer model to be valid at the chord, not merely flat."""
    pytest.importorskip("ngsolve")
    from netgen.occ import MoveTo

    from radia.sibc_corner_patch import solve_cut_patch

    # 0.2 mm thick at 150 kHz: about one skin depth, not the six required.
    half = 0.1e-3
    corners = np.array([[-2e-3, -half], [2e-3, -half],
                        [2e-3, half], [-2e-3, half]])
    dense = np.concatenate([
        np.linspace(corners[i], corners[(i + 1) % 4], 64, endpoint=False)
        for i in range(4)])
    ds = np.linalg.norm(dense - np.roll(dense, 1, axis=0), axis=1)
    face = MoveTo(-2e-3, -half).Rectangle(4e-3, 2 * half).Face()

    with pytest.raises(ValueError, match="skin depths thick"):
        solve_cut_patch(
            face, x_cut=0.0, panel_xy=dense, panel_ds=ds,
            panel_current=np.full(len(dense), 1.0 / len(dense), dtype=complex),
            axial_E_field=1.0 + 0j, surface_impedance=1e-3 + 1e-3j,
            omega=2.0 * math.pi * 150_000.0, sigma=5.8e7,
            maxh_conductor=2e-5)


def test_surface_potential_inverts_the_panel_equation():
    """``A`` at the surface comes from the outer solve's own equation.

    ``Zs K + j omega A = E0`` holds panel by panel, so the surface potential
    is already in the outer solution and does not have to be recovered by
    re-integrating a singular kernel over the sheet the point sits on.
    """
    from radia.sibc_corner_patch import surface_potential

    rng = np.random.default_rng(7)
    ds = rng.uniform(1e-5, 3e-5, 32)
    current = rng.normal(size=32) + 1j * rng.normal(size=32)
    zs, e0, omega = 3.1e-4 + 3.1e-4j, 0.25 - 0.1j, 2 * math.pi * 150_000.0

    a = surface_potential(current, ds, surface_impedance=zs,
                          axial_E_field=e0, omega=omega)
    residual = zs * (current / ds) + 1j * omega * a - e0
    assert np.max(np.abs(residual)) < 1e-12 * abs(e0)
