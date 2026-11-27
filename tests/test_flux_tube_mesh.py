"""Surrogate contracts for the experimental conforming hex flux-tube mesh.

Guards: ring equal-flux cells (counts, shared sector breaks, cell flux), the
spherical cap chart and its two-crossing chords, the Bernstein Jacobian bound
(also far from the origin), and on one coarse sphere
case with a 2:1 rim: face conformity, the 2:1 seam (every fine face shared
once, opposite orientation), positive Jacobian bound inside every hex,
section flux against the analytic tube flux, and native Gmsh readback of the
written MSH. A straight-wire case guards the planar balanced route.
Assumes closed forms: a uniformly magnetized sphere (mu0 M = 1 T along z) has
B = 2/3 T inside, so a disk or cap of radius R inside carries 2/3 pi R^2 Wb;
a straight wire carries mu0 I h ln(r2 / r1) / (2 pi) through a radial
rectangle. The flux bounds are coarse-mesh regression guards (bilinear
section quads, 4 segments per cell side), not accuracy criteria.
"""

import re
import shutil
import subprocess

import numpy as np
import pytest

from radia.flux_seeding import (annular_patch, balanced_equal_flux_seeds, planar_patch,
                                ring_equal_flux_cells, spherical_cap_patch,
                                surface_crossings)
from radia.flux_tube_mesh import (STATUS, build_flux_tube_mesh, face_census, flux_tube_report,
                                  hex_corner_jacobians, hex_jacobian_lower_bound,
                                  write_flux_tube_msh)

MU0 = 4e-7 * np.pi
A = 0.05
UNIT_CUBE = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                      [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], float)


def _sphere_field(radius, b_r):
    m = 4 / 3 * np.pi * radius ** 3 * b_r / MU0

    def field(p):
        p = np.asarray(p, dtype=float)
        shape = p.shape
        p = p.reshape(-1, 3)
        r = np.linalg.norm(p, axis=1)
        out = np.tile([0.0, 0.0, 2 / 3 * b_r], (len(p), 1))
        o = r >= radius
        rh = p[o] / r[o][:, None]
        out[o] = MU0 * m / (4 * np.pi) * (3 * rh[:, 2:3] * rh - [0, 0, 1]) / r[o][:, None] ** 3
        return out.reshape(shape)
    return field


def _wire_field(current):
    def field(p):
        p = np.asarray(p, dtype=float)
        shape = p.shape
        p = p.reshape(-1, 3)
        k = MU0 * current / (2 * np.pi * (p[:, 0] ** 2 + p[:, 1] ** 2))
        return np.column_stack((-k * p[:, 1], k * p[:, 0], 0 * k)).reshape(shape)
    return field


def _meridian_angle(p):
    """Angle about the equator circle (rho, z) = (A, 0), z stretched near it."""
    p = np.asarray(p, float)
    z = p[..., 2]
    z = z * np.abs(z) / (np.abs(z) + 2 * A)
    return np.arctan2(z, np.hypot(p[..., 0], p[..., 1]) - A)


def test_ring_cells_carry_equal_flux_and_share_sector_breaks():
    patch = annular_patch([0, 0, 0], [0, 0, 1], 0.04)
    field = lambda p: np.broadcast_to([0.0, 0.0, 0.3], np.shape(p))  # noqa: E731
    cells = ring_equal_flux_cells(patch.flux_density(field), patch.metric, n_lines=40,
                                  rim_factor=2, transition=True)
    nc, nx, ny = cells.n_centre, cells.n_rings, cells.n_sectors
    assert nc + ny * (nx - 1 + 2) == 40, (nc, nx, ny)
    assert ny % 2 == 0 and ny % nc == 0, ny
    assert cells.transition_cells == tuple(range(nc + (nx - 2) * ny, nc + (nx - 1) * ny))
    c = cells.cells
    rim = c[nc + (nx - 1) * ny:]
    ring_breaks = np.unique(c[nc:nc + ny, 2:].ravel())
    assert np.all(np.isin(ring_breaks, np.unique(rim[:, 2:].ravel()))), "rim must contain ring breaks"
    r0, r1 = 0.04 * c[:, 0], 0.04 * c[:, 1]
    flux = 0.3 * (r1 ** 2 - r0 ** 2) / 2 * 2 * np.pi * (c[:, 3] - c[:, 2])
    target = 0.3 * np.pi * 0.04 ** 2 / 40
    assert np.allclose(flux, target, rtol=2e-3), (flux.min() / target, flux.max() / target)
    assert cells.target_cell_flux_Wb == pytest.approx(target, rel=1e-6)


def test_ring_cells_reject_impossible_layouts():
    patch = annular_patch([0, 0, 0], [0, 0, 1], 0.04)
    g = patch.flux_density(lambda p: np.broadcast_to([0.0, 0.0, 1.0], np.shape(p)))
    with pytest.raises(ValueError, match="no ring layout"):
        ring_equal_flux_cells(g, patch.metric, n_lines=7)
    with pytest.raises(ValueError, match="rim_factor = 2"):
        ring_equal_flux_cells(g, patch.metric, n_lines=40, rim_factor=3, transition=True)


def test_spherical_cap_chart_area_and_locate():
    cap = spherical_cap_patch([0.1, 0, 0], [0, 0, 1], 0.03, np.pi / 3)
    t, w = np.polynomial.legendre.leggauss(12)
    t, w = 0.5 * (t + 1), 0.5 * w
    XI, ET = np.meshgrid(t, t, indexing="ij")
    g = cap.flux_density(lambda p: np.broadcast_to([0.0, 0.0, 2.0], np.shape(p)))
    flux = float(np.sum(np.outer(w, w) * g(XI, ET)))
    assert flux == pytest.approx(2.0 * np.pi * (0.03 * np.sin(np.pi / 3)) ** 2, rel=1e-10)
    s, xi, eta = cap.locate(cap.chart(XI, ET))
    assert np.max(np.abs(s)) < 1e-15 and np.allclose(xi, XI) and np.allclose(eta, ET)


def test_a_chord_through_a_sphere_cap_crosses_it_twice():
    """Guards: the cap splits every segment at its point nearest the centre.

    Both ends of the segment x = -2 -> 2 (z = 0.5) lie outside the unit sphere, so
    endpoint signs alone see no crossing; the analytic roots x = -+sqrt(0.75) lie on the
    upper hemisphere, entering (-1) then leaving (+1). Splitting the chord at x = 0
    by hand must give the same answer.
    """
    cap = spherical_cap_patch([0, 0, 0], [0, 0, 1], 1.0, np.pi / 2)
    chord = np.array([[-3.0, 0, 0.5], [-2.0, 0, 0.5], [2.0, 0, 0.5]])
    split = np.array([[-3.0, 0, 0.5], [-2.0, 0, 0.5], [0.0, 0, 0.5], [2.0, 0, 0.5]])
    assert surface_crossings(chord, cap, seed_vertex=0, closed=False) == [-1, 1]
    assert surface_crossings(split, cap, seed_vertex=0, closed=False) == [-1, 1]
    below = chord - [0, 0, 1.0]                     # z = -0.5: lower hemisphere only
    assert surface_crossings(below, cap, seed_vertex=0, closed=False) == []


def _sampled_jacobians(X, points):
    """Direct det(dx/du) of one Gmsh-ordered trilinear hex at reference points."""
    out = []
    for q in np.asarray(points, float):
        f = np.where(UNIT_CUBE > 0, q, 1 - q)
        d = np.stack([(2 * UNIT_CUBE[:, a] - 1) * np.prod(np.delete(f, a, axis=1), axis=1)
                      for a in range(3)], 1)                      # (8 nodes, 3 ref axes)
        out.append(np.linalg.det((X - X[0]).T @ d))
    return np.array(out)


def test_jacobian_bound_sees_an_interior_inversion_far_from_the_origin():
    """Guards: the bound recentres each hex and keeps a rounding margin.

    All corner determinants are positive, but det < 0 near reference (0.2, 0, 0); far
    from the origin (1e15 m, spacing 1/8 representable) a bound formed from absolute
    coordinates came out positive. Oracle: direct sampling of the derivative.
    """
    X = np.array([[.75, .125, .375], [.875, -.5, .75], [1, 1, .25], [-.125, 1.25, .375],
                  [0, 0, 1.25], [1.25, 0, 1], [1.375, 1.125, 1.125], [.25, 1.25, 1.375]])
    hexes = np.arange(8)[None]
    probe = [[0.2, 0.0, 0.0], [0.2, 1e-4, 1e-4]]
    for shift in (0.0, 1e15):
        Y = X + shift
        assert np.all(hex_corner_jacobians(Y, hexes) > 0)
        direct = _sampled_jacobians(Y, probe)
        bound = hex_jacobian_lower_bound(Y, hexes)[0]
        assert direct.max() < 0, f"shift {shift}: sampled det {direct}"
        assert bound <= direct.min(), f"shift {shift}: bound {bound} above sampled {direct}"
    # Recentring keeps the bound useful far away: the unit cube still certifies det = 1.
    assert hex_jacobian_lower_bound(UNIT_CUBE + 1e15, hexes)[0] == pytest.approx(1.0, rel=1e-9)


def test_jacobian_bound_never_certifies_a_singular_hex():
    """Guards: the rounding margin. Singular affine hexes have det = 0 exactly everywhere,
    so the computed Bernstein minimum is rounding noise and sometimes slightly positive;
    with the margin no such hex may be reported positive."""
    rng = np.random.default_rng(5)
    n = 2000
    A = rng.standard_normal((n, 3, 3))
    A[:, :, 2] = A[:, :, 0] * np.sqrt(2) - A[:, :, 1] / 3
    X = np.einsum("vk,nik->nvi", UNIT_CUBE, A) + rng.standard_normal((n, 1, 3))
    bound = hex_jacobian_lower_bound(X.reshape(-1, 3), np.arange(8 * n).reshape(n, 8))
    assert np.count_nonzero(bound > 0) == 0, f"max bound {bound.max()}"

def test_jacobian_bound_is_exact_for_affine_and_below_sampled_minimum():
    M = np.array([[2.0, 0.3, 0.0], [0.0, 1.0, 0.2], [0.1, 0.0, 0.5]])
    hexes = np.arange(8)[None]
    assert hex_jacobian_lower_bound(UNIT_CUBE @ M.T, hexes)[0] == pytest.approx(np.linalg.det(M))
    rng = np.random.default_rng(3)
    t = np.linspace(0, 1, 11)
    U, V, W = (x.ravel() for x in np.meshgrid(t, t, t, indexing="ij"))
    for _ in range(20):
        X = UNIT_CUBE + 0.3 * rng.standard_normal((8, 3))
        bound = hex_jacobian_lower_bound(X, hexes)[0]
        f = [[(U if r else 1 - U), (V if s else 1 - V), (W if q else 1 - W)] for r, s, q in UNIT_CUBE]
        d = [[(1 if r else -1) * f[n][1] * f[n][2], (1 if s else -1) * f[n][0] * f[n][2],
              (1 if q else -1) * f[n][0] * f[n][1]] for n, (r, s, q) in enumerate(UNIT_CUBE)]
        J = np.einsum("nak,ni->kia", np.array(d), X)
        assert bound <= np.linalg.det(J).min() + 1e-12


@pytest.fixture(scope="module")
def sphere_mesh():
    field = _sphere_field(A, 1.0)
    patch = annular_patch([0, 0, 0], [0, 0, 1], 0.999 * A)
    seeds = ring_equal_flux_cells(patch.flux_density(field), patch.metric, n_lines=16,
                                  rim_factor=2, transition=True)
    mesh = build_flux_tube_mesh(
        field, patch, seeds, angle_fn=_meridian_angle, n_sections=8, step_m=0.002,
        max_steps=3000, extent_m={"min": [-1.0] * 3, "max": [1.0] * 3}, cells_per_side=4,
        periodic_eta=True, collapsed_xi0=True, rim_grading=(1 / 0.999, 4),
        interface=([0.0, 0.0, 0.0], A), radial_centre=[0.0, 0.0, 0.0],
        level_target_radius=0.9 * 12.0 * A)
    return mesh, seeds


def test_sphere_tubes_conform_and_stay_valid(sphere_mesh):
    mesh, _ = sphere_mesh
    rep = flux_tube_report(mesh)
    assert rep["status"] == STATUS == "experimental"
    assert rep["untraced_tubes"] == 0, rep
    assert rep["nonpositive_jacobian_bound_hex"] == 0, rep["min_jacobian_bound_m3"]
    faces = rep["faces"]
    assert faces["faces_used_more_than_twice"] == 0, faces
    assert faces["interior_faces_same_orientation"] == 0, faces
    assert faces["boundary_faces_on_another_tube"] == 0, faces
    assert mesh.info["closed_walls_not_one_turn"] == 0, mesh.info


def test_two_to_one_seam_shares_every_fine_face_once(sphere_mesh):
    mesh, seeds = sphere_mesh
    trans = set(seeds.transition_cells)
    nc, nx, ny = seeds.n_centre, seeds.n_rings, seeds.n_sectors
    rim = set(range(nc + (nx - 1) * ny, len(seeds.cells)))
    census = face_census(mesh.hexes, mesh.tube_of)
    a, b = census["shared_face_pairs"]
    ta, tb = mesh.tube_of[a], mesh.tube_of[b]
    seam = (np.isin(ta, list(trans)) & np.isin(tb, list(rim))) | (np.isin(tb, list(trans)) & np.isin(ta, list(rim)))
    expected = 0
    for s in mesh.tubes:
        if s["cell"] in trans:
            k0, k1 = s["k_range"]
            layers = (k1 - k0 + 1) if s["closed"] else (k1 - k0)
            expected += s["fine_outer_segments"] * layers
            assert s["fine_outer_segments"] == 2 * s["grid_nx_ny"][1], s
    assert expected > 0
    assert int(np.sum(seam)) == expected, (int(np.sum(seam)), expected)
    assert census["interior_faces_same_orientation"] == 0


def test_sphere_section_flux_matches_analytic_tube_flux(sphere_mesh):
    mesh, seeds = sphere_mesh
    analytic = 2 / 3 * np.pi * (0.999 * A) ** 2 / 16
    assert seeds.target_cell_flux_Wb == pytest.approx(analytic, rel=1e-3)
    flux = np.concatenate([s["section_flux_Wb"] for s in mesh.tubes]) / analytic
    assert 0.95 < flux.min() and flux.max() < 1.05, (flux.min(), flux.max())


def _readback(path, tmp_path):
    """(node count, element type codes, node coordinates) as Gmsh itself reads the file."""
    try:
        import gmsh
    except ImportError:
        gmsh = None
    if gmsh is not None:
        gmsh.initialize()
        try:
            gmsh.open(str(path))
            tags, xyz, _ = gmsh.model.mesh.getNodes()
            types = list(gmsh.model.mesh.getElements(3)[0])
            n_el = sum(len(e) for e in gmsh.model.mesh.getElements(3)[1])
            return len(tags), types, n_el, np.asarray(xyz).reshape(-1, 3)[np.argsort(tags)]
        finally:
            gmsh.finalize()
    exe = shutil.which("gmsh")
    if exe is None:
        pytest.skip("Gmsh is not available (neither the Python API nor gmsh on PATH)")
    out = tmp_path / "readback.msh"
    subprocess.run([exe, str(path), "-save", "-format", "msh22", "-o", str(out)],
                   check=True, capture_output=True, timeout=120)
    text = out.read_text(encoding="utf-8")
    nodes = re.search(r"\$Nodes\n(\d+)\n(.*?)\$EndNodes", text, re.S)
    rows = np.array([list(map(float, ln.split())) for ln in nodes.group(2).splitlines()])
    els = re.search(r"\$Elements\n(\d+)\n(.*?)\$EndElements", text, re.S)
    types = sorted({int(ln.split()[1]) for ln in els.group(2).splitlines()})
    return int(nodes.group(1)), types, int(els.group(1)), rows[np.argsort(rows[:, 0]), 1:]


def test_msh_reads_back_in_gmsh(sphere_mesh, tmp_path):
    mesh, _ = sphere_mesh
    path = tmp_path / "tubes.msh"
    write_flux_tube_msh(path, mesh, length_scale=1e-3)
    n_nodes, types, n_el, xyz = _readback(path, tmp_path)
    assert n_nodes == len(mesh.nodes_m) and n_el == len(mesh.hexes), (n_nodes, n_el)
    assert list(types) == [5], types
    assert np.allclose(xyz, mesh.nodes_m / 1e-3, rtol=0, atol=1e-9)


def test_wire_tubes_on_balanced_cells():
    current, r1, r2, h = 10.0, 0.01, 0.05, 0.04
    field = _wire_field(current)
    patch = planar_patch([r1, 0, -h / 2], [r2 - r1, 0, 0], [0, 0, h])
    nodes = np.linspace(0, 1, 33)
    seeds = balanced_equal_flux_seeds(patch.flux_density(field), xi_nodes=nodes, eta_nodes=nodes,
                                      n_lines=12, chart=patch.chart, metric=patch.metric)
    mesh = build_flux_tube_mesh(
        field, patch, seeds, angle_fn=lambda p: np.arctan2(np.asarray(p)[..., 1], np.asarray(p)[..., 0]),
        n_sections=12, step_m=0.001, max_steps=1000, cells_per_side=3)
    rep = flux_tube_report(mesh)
    assert rep["untraced_tubes"] == 0 and rep["nonpositive_jacobian_bound_hex"] == 0, rep
    assert rep["faces"]["faces_used_more_than_twice"] == 0
    assert rep["faces"]["interior_faces_same_orientation"] == 0
    assert rep["faces"]["boundary_faces_on_another_tube"] == 0
    analytic = MU0 * current * h * np.log(r2 / r1) / (2 * np.pi) / 12
    assert abs(seeds.target_cell_flux_Wb) == pytest.approx(analytic, rel=1e-3)
    sign = np.sign(seeds.target_cell_flux_Wb)       # the patch normal is -y here
    flux = sign * np.concatenate([s["section_flux_Wb"] for s in mesh.tubes]) / analytic
    assert 0.95 < flux.min() and flux.max() < 1.05, (flux.min(), flux.max())
    assert np.all(hex_corner_jacobians(mesh.nodes_m, mesh.hexes) > 0)


def test_ring_centre_cells_get_one_even_eta_count():
    """Guards: collapsed centre cells whose wall nodes give an odd or unequal eta count are
    brought to one even count (halving their longest segments) instead of raising.
    Assumes: an iron ring with a toroidal winding, B = R/rho e_phi (R = 0.1 m, minor radius
    b = 0.03 m), so the flux through a half-plane section is R 2 pi (R - sqrt(R^2 - b^2))."""
    R, b = 0.1, 0.03

    def field(p):
        p = np.asarray(p, float)
        s = R / (p[..., 0] ** 2 + p[..., 1] ** 2)
        return np.stack([-p[..., 1] * s, p[..., 0] * s, 0 * s], axis=-1)

    patch = annular_patch([R, 0, 0], [0, 1, 0], b, r_inner_m=0.0)
    seeds = ring_equal_flux_cells(patch.flux_density(field), patch.metric, n_lines=8)
    mesh = build_flux_tube_mesh(
        field, patch, seeds, angle_fn=lambda p: np.arctan2(np.asarray(p)[..., 1], np.asarray(p)[..., 0]),
        n_sections=4, step_m=0.002, max_steps=3000, extent_m={"min": [-1.0] * 3, "max": [1.0] * 3},
        cells_per_side=4, periodic_eta=True, collapsed_xi0=True)
    rep = flux_tube_report(mesh)
    assert rep["untraced_tubes"] == 0 and rep["nonpositive_jacobian_bound_hex"] == 0, rep
    assert rep["faces"]["faces_used_more_than_twice"] == 0, rep["faces"]
    assert rep["faces"]["interior_faces_same_orientation"] == 0, rep["faces"]
    assert rep["faces"]["boundary_faces_on_another_tube"] == 0, rep["faces"]
    assert mesh.info["closed_walls_not_one_turn"] == 0, mesh.info
    ny = {s["grid_nx_ny"][1] for s in mesh.tubes if s["cell"] < seeds.n_centre}
    assert len(ny) == 1 and next(iter(ny)) % 2 == 0, ny
    exact = R * 2 * np.pi * (R - np.sqrt(R * R - b * b))
    total = np.sum([s["section_flux_Wb"] for s in mesh.tubes], axis=0) / exact
    assert np.all(np.abs(total - 1) < 0.03), total
