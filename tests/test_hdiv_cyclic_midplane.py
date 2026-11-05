"""A connected sector with the midplane mirror reproduces the closed ring, on the CPU route.

Guards: `image="-z"` together with `image_cyclic` and the azimuthal periodic trace.
Assumes: an annular HEX yoke, a field along +z, a mesh that is the mirror of its upper half.
Not independent validation: the ring is the same route with every image removed.
"""
import math

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia import vim
from radia.vim import _image

SEAM = ("sector_low", "sector_high")
MU_R = 10.0
FOLD = 8
LIMIT = 1e-6
PROBES = ((0.024, 0.0, 0.006), (0.04, 0.0, 0.006), (0.06, 0.0, 0.006), (0.04, 0.01, 0.009))


def _yoke(path, cells, fold, levels, *, ring=False, inner=0.03, outer=0.05):
    """Annular HEX yoke with its height split at `levels`: the closed ring, or one sector with identified seam faces."""
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, IdentificationType, Mesh, MeshPoint, Pnt
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "yoke")
    skin = mesh.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    mesh.SetBCName(0, "skin")
    if not ring:
        low, high = mesh.Add(FaceDescriptor(surfnr=2, domin=1, domout=0, bc=2)), mesh.Add(FaceDescriptor(surfnr=3, domin=1, domout=0, bc=3))
        mesh.SetBCName(1, SEAM[0])
        mesh.SetBCName(2, SEAM[1])
    total = cells*fold if ring else cells
    stations = total if ring else total + 1
    step = 2.0*math.pi/(cells*fold)
    point = {}
    for it, ir, iz in np.ndindex(stations, 2, len(levels)):
        radius = (inner, outer)[ir]
        point[it, ir, iz] = mesh.Add(MeshPoint(Pnt(radius*math.cos(it*step), radius*math.sin(it*step), levels[iz])))
    for it, iz in np.ndindex(total, len(levels) - 1):
        following = (it + 1) % stations
        corners = [point[it, 0, iz], point[it, 1, iz], point[following, 1, iz], point[following, 0, iz],
                   point[it, 0, iz + 1], point[it, 1, iz + 1], point[following, 1, iz + 1], point[following, 0, iz + 1]]
        mesh.Add(Element3D(1, corners))
        quads = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (2, 3, 7, 6), (0, 4, 7, 3), (1, 2, 6, 5)]
        for index, quad in enumerate(quads):
            if index == 0 and iz != 0 or index == 1 and iz != len(levels) - 2:
                continue
            if index in (2, 3) and not ring and (index == 2 and it == 0 or index == 3 and it == total - 1):
                mesh.Add(Element2D(low if index == 2 else high, [corners[k] for k in quad]))
            elif index < 2 or index > 3:
                mesh.Add(Element2D(skin, [corners[k] for k in quad]))
    if not ring:
        for ir, iz in np.ndindex(2, len(levels)):
            mesh.AddPointIdentification(point[0, ir, iz], point[total, ir, iz], identnr=1, type=IdentificationType.PERIODIC)
    mesh.Save(str(path))
    return ng.Mesh(str(path))


def _solve(mesh, order, **images):
    with ng.TaskManager():
        return vim.Solve(mesh, mu_r=MU_R, H_ext=ng.CoefficientFunction((0, 0, 1e5)), order=order, tol=1e-10, **images)


@pytest.fixture(scope="module")
def solved(tmp_path_factory):
    folder = tmp_path_factory.mktemp("cyclic_midplane")
    results = {}
    for order in (1, 2):
        ring = _yoke(folder/f"ring{order}.vol", 2, FOLD, (-0.004, 0.0, 0.004), ring=True)
        half = _yoke(folder/f"half{order}.vol", 2, FOLD, (0.0, 0.004))
        full = _yoke(folder/f"full{order}.vol", 2, FOLD, (-0.004, 0.0, 0.004))
        results[order] = dict(ring=_solve(ring, order), half=_solve(half, order, image="-z", image_cyclic=FOLD, cyclic_periodic_boundaries=SEAM),
                              full=_solve(full, order, image_cyclic=FOLD, cyclic_periodic_boundaries=SEAM),
                              lone=_solve(half, order))
    return results


@pytest.mark.parametrize("order", (1, 2))
def test_the_half_sector_reproduces_the_ring_and_the_full_height_sector(solved, order):
    ring, half, full, lone = (solved[order][name] for name in ("ring", "half", "full", "lone"))
    reference = ring["M_avg"][2]
    for name, value in (("half", half["M_avg"][2]), ("full", full["M_avg"][2])):
        assert abs(value - reference) <= LIMIT*abs(reference), f"order={order} {name} Mz={value!r} ring Mz={reference!r} limit={LIMIT:.0e}"
    assert np.max(np.abs(half["M_avg"][:2])) <= LIMIT*abs(reference), f"order={order} in-plane mean {half['M_avg']!r} of the full domain"
    assert abs(lone["M_avg"][2] - reference) > 1e3*LIMIT*abs(reference), f"a half sector without images passes: {lone['M_avg']!r} ring {reference!r}"


def test_a_plus_midplane_or_another_mirror_with_rotations_refuses():
    assert _image.cyclic_image_terms([(2, -1)], 2)[0] == [0, 4, 4]
    for planes in ([(2, 1)], [(0, -1)], [(2, -1), (0, -1)]):
        with pytest.raises(NotImplementedError, match="midplane"):
            _image.cyclic_image_terms(planes, 4)


@pytest.mark.parametrize("fold,alternating", ((4, False), (6, True)))
def test_the_term_list_is_closed_under_inversion_with_matching_signs(fold, alternating):
    masks, signs, angles = _image.cyclic_image_terms([(2, -1)], fold, alternating)
    assert len(masks) == 2*fold - 1, (len(masks), fold)
    assert angles[0] != 0.0, f"the seam pairing reads the first angle: {angles!r}"
    rows = {(mask, round(math.remainder(angle, 2.0*math.pi), 12)): sign for mask, sign, angle in zip(masks, signs, angles)}
    for (mask, angle), sign in rows.items():
        partner = (mask, round(math.remainder(-angle, 2.0*math.pi), 12))
        partner = partner if partner in rows else (mask, -partner[1])
        assert rows[partner] == sign, (mask, angle, sign, rows[partner])
