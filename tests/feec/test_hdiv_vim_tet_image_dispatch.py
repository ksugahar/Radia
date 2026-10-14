"""Mirror-image charge Gram on flat and curved TET meshes: the host-block / far dispatch.

A quarter model (x > 0, y > 0) with ``image="-x-y"`` -- the quadrupole parity, potential
``-g x y`` odd in both x and y -- must reproduce its full model.  For the flat box the full mesh is
the exact mirrored union of the quarter mesh, so the reduced+image solve and the full solve differ
only by ACA truncation, Krylov tolerance, and summation order.  For the curved sphere the two
meshes are independent, so the agreement is at discretization level.

The folded Gram must also stay what the Coulomb energy makes it: symmetric, positive definite, and
zero on the cut-face charges that lie on an antisymmetric mirror plane (their image annihilates
them).  These checks run on the sigma-normalized dense Gram assembled from the symmetric H-matrix
apply, with the raw O(n^2) quadratic form confirming the sign without the H-matrix in the loop.

The dispatch counters prove the far / host-block rules were actually used: before them every
mirror image on a curved mesh was a per-entry scalar curved Duffy (validation_test/feec/
validate_hdiv_vim_tet_image_dispatch.py holds the legacy A/B and the timing).
"""

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
from netgen.meshing import (  # noqa: E402
    Element2D, Element3D, FaceDescriptor, Mesh as NetgenMesh, MeshPoint, Pnt as MeshPnt)
from netgen.occ import Box, OCCGeometry, Pnt, Sphere  # noqa: E402

from radia import vim  # noqa: E402
from radia.vim._image import image_group, parse_image_string  # noqa: E402

MU_R = 100.0
GRADIENT = 1.0e4
IMAGE = "-x-y"
OBSERVATION_POINTS = np.array([
    [1.5, 0.5, 0.2], [-1.5, 0.5, 0.2], [0.5, -1.5, 0.2], [-0.5, -1.5, -0.2],
    [1.2, 1.2, 0.5], [-1.2, -1.2, -0.5], [0.3, 0.3, 1.6], [-0.3, 0.3, -1.6],
])


def _quadrupole():
    return ng.CoefficientFunction((GRADIENT * ng.y, GRADIENT * ng.x, 0.0))


def _quarter_box(maxh):
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Box(Pnt(0, 0, -1), Pnt(1, 1, 1))).GenerateMesh(maxh=maxh))


def _mirrored_union(quarter, axes=(0, 1), plane_tol=1.0e-12):
    source = quarter.ngmesh
    coordinates = [tuple(float(v) for v in point.p) for point in source.Points()]
    full = NetgenMesh(dim=3)
    full.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    subsets = [()] + [tuple(axes[k] for k in range(len(axes)) if mask & (1 << k))
                      for mask in range(1, 1 << len(axes))]
    index, new_id = {}, {}
    for subset in subsets:
        for old, xyz in enumerate(coordinates, start=1):
            reflected = list(xyz)
            for axis in subset:
                reflected[axis] = -reflected[axis]
            key = tuple(round(v, 12) + 0.0 for v in reflected)
            if key not in index:
                index[key] = full.Add(MeshPoint(MeshPnt(*reflected)))
            new_id[(subset, old)] = index[key]
    for subset in subsets:
        flip = len(subset) % 2 == 1
        for element in source.Elements3D():
            vertices = [new_id[(subset, v.nr)] for v in element.vertices]
            if flip:
                vertices[0], vertices[1] = vertices[1], vertices[0]
            full.Add(Element3D(1, vertices))
        for element in source.Elements2D():
            if any(all(abs(coordinates[v.nr - 1][axis]) <= plane_tol for v in element.vertices)
                   for axis in axes):
                continue
            vertices = [new_id[(subset, v.nr)] for v in element.vertices]
            if flip:
                vertices[0], vertices[1] = vertices[1], vertices[0]
            full.Add(Element2D(1, vertices))
    full.SetMaterial(1, "iron")
    return ng.Mesh(full)


def _quarter_sphere(maxh):
    shape = Sphere(Pnt(0, 0, 0), 1.0) * Box(Pnt(0, 0, -2), Pnt(2, 2, 2))
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(shape).GenerateMesh(maxh=maxh))
        mesh.Curve(2)
    return mesh


def _full_sphere(maxh):
    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=maxh))
        mesh.Curve(2)
    return mesh


def _field(mesh, order, curved, image):
    with ng.TaskManager():
        result = vim.Solve(
            mesh, mu_r=MU_R, H_ext=_quadrupole(), order=order,
            curve_order=2 if curved else None, gram_eps=1e-12, tol=1e-12, image=image)
        return np.asarray(vim.FieldFromSolution(result, OBSERVATION_POINTS, algorithm="direct"))


def _relative_error(reference, candidate):
    scale = np.maximum(np.linalg.norm(reference, axis=1),
                       GRADIENT * np.linalg.norm(OBSERVATION_POINTS[:, :2], axis=1))
    return float(np.max(np.linalg.norm(candidate - reference, axis=1) / scale))


def _folded_gram(mesh, order, curved):
    masks, signs = [], []
    for axes, sign in image_group(parse_image_string(IMAGE)):
        masks.append(int(sum(1 << axis for axis in axes)))
        signs.append(float(sign))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        B, G, _ = vim.ChargeGram(
            fes, eps=1e-12, curve_order=2 if curved else None,
            image_masks=masks, image_signs=signs)
    return fes, B.tocsr(), G


def _plane_charge_rows(mesh, fes, B, axes=(0, 1), plane_tol=1.0e-12):
    plane_dofs = set()
    for element in mesh.Elements(ng.BND):
        coordinates = np.asarray([mesh.vertices[v.nr].point for v in element.vertices], float)
        if any(np.max(np.abs(coordinates[:, axis])) <= plane_tol for axis in axes):
            plane_dofs.update(int(d) for d in fes.GetDofNrs(element) if int(d) >= 0)
    rows = []
    for row in range(B.shape[0]):
        columns = B.indices[B.indptr[row]:B.indptr[row + 1]]
        if len(columns) and all(int(c) in plane_dofs for c in columns):
            rows.append(row)
    return rows


@pytest.mark.parametrize("order", (1, 2))
def test_flat_box_quarter_with_image_reproduces_mirrored_full_model(order):
    quarter = _quarter_box(0.55)
    full = _mirrored_union(quarter)
    assert full.ne == 4 * quarter.ne
    reduced_field = _field(quarter, order, False, IMAGE)
    full_field = _field(full, order, False, None)
    assert _relative_error(full_field, reduced_field) < 1.0e-8


@pytest.mark.parametrize("order", (1, 2))
def test_curved_sphere_quarter_with_image_matches_full_sphere(order):
    reduced_field = _field(_quarter_sphere(0.55), order, True, IMAGE)
    full_field = _field(_full_sphere(0.55), order, True, None)
    assert _relative_error(full_field, reduced_field) < 3.0e-2


@pytest.mark.parametrize("geometry,order", [("flat_box", 1), ("flat_box", 2),
                                            ("curved_sphere", 1), ("curved_sphere", 2)])
def test_folded_gram_is_symmetric_positive_and_annihilates_plane_charges(geometry, order, monkeypatch):
    monkeypatch.setenv("RADIA_HDIV_BLOCK_CACHE_STATS", "1")
    curved = geometry == "curved_sphere"
    mesh = _quarter_sphere(0.6) if curved else _quarter_box(0.6)
    fes, B, G = _folded_gram(mesh, order, curved)
    n = int(G.ndof())
    with ng.TaskManager():
        diagonal = np.array([G.entry(p, p) for p in range(n)])
        rng = np.random.default_rng(20260905)
        pairs = rng.integers(0, n, size=(600, 2))
        forward = np.array([G.entry(int(a), int(b)) for a, b in pairs])
        backward = np.array([G.entry(int(b), int(a)) for a, b in pairs])
        sigma = np.asarray(G.charge_sigma(), float)
        dense = np.empty((n, n))
        unit = np.zeros(n)
        for j in range(n):
            unit[:] = 0.0
            unit[j] = 1.0
            dense[:, j] = np.asarray(G.matvec_sym(unit), float)
        normalized = dense / sigma[:, None] / sigma[None, :]
        eigenvalues, vectors = np.linalg.eigh(0.5 * (normalized + normalized.T))
        raw_form = float(G.raw_symmetric_quadratic_form(vectors[:, 0] / sigma))
        stats = dict(G.stats())
    scale = float(np.max(np.abs(diagonal)))
    assert np.max(np.abs(forward - backward)) / scale < 1.0e-12
    plane_rows = _plane_charge_rows(mesh, fes, B)
    assert plane_rows, "the quarter mesh must carry cut-face charges on the mirror planes"
    other = np.setdiff1d(np.arange(n), np.array(plane_rows))
    assert np.max(np.abs(diagonal[plane_rows])) / np.median(diagonal[other]) < 1.0e-10
    assert np.all(sigma[plane_rows] == 1.0)
    assert np.max(np.abs(normalized - normalized.T)) < 1.0e-12
    assert eigenvalues[0] > -1.0e-9
    assert raw_form > -1.0e-9
    assert stats["ho_image_far_entries"] > 0.0
    assert stats["ho_image_far_enabled"] == 1.0
    assert stats["ho_image_block_enabled"] == 1.0
    if curved or order >= 2:
        # A host block exists: curved always, flat only with the analytic BDM2 block.
        assert stats["ho_image_block_entries"] > 0.0
        assert stats["ho_image_scalar_entries"] == 0.0
    else:
        # Flat BDM1 keeps the per-entry analytic fold for near image pairs by design;
        # only the far rule changes there.
        assert stats["ho_image_block_entries"] == 0.0
        assert stats["ho_image_scalar_entries"] > 0.0
