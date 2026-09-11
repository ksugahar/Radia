"""Which faces add_kelvin_cubit treats as the Kelvin interface, without Cubit.

A fake ``cubit`` object answers the few queries the selection and pairing
helpers make.  The regression: the largest face of EVERY air volume used to go
into kelvin_int, so the faces of gap slabs inside the air did too, and the
exterior mesh was paired by ``min(len(air), len(kelvin))``, which silently
dropped any face without a partner.
"""
import importlib.util
import math
from pathlib import Path

import pytest

path = (Path(__file__).resolve().parents[1] / "packages" / "cubit-mesh-export"
        / "src" / "cubit_mesh_export" / "cubit_helpers" / "add_kelvin.py")
spec = importlib.util.spec_from_file_location("add_kelvin_under_test", path)
add_kelvin = importlib.util.module_from_spec(spec)
spec.loader.exec_module(add_kelvin)

R = 0.22
OFFSET = (3 * R, 0.0, 0.0)
ORIGIN = (0.0, 0.0, 0.0)


class SphereFace:
    """A spherical cap {|x - c| = radius, z_low <= (x - c)_z <= z_high}."""

    def __init__(self, centre, radius=R, z_low=None, z_high=None):
        self.centre, self.radius = centre, radius
        self.z_low, self.z_high = z_low, z_high

    def closest_point_trimmed(self, point):
        d = [point[i] - self.centre[i] for i in range(3)]
        norm = math.sqrt(sum(v * v for v in d))
        if norm == 0.0:
            # every point is closest to the centre; Cubit returns one on the
            # face (measured: distance R), so pick one on the equator
            d, norm = [1.0, 0.0, 0.0], 1.0
        u = [v / norm * self.radius for v in d]
        clamp = None
        if self.z_low is not None and u[2] < self.z_low:
            clamp = self.z_low
        if self.z_high is not None and u[2] > self.z_high:
            clamp = self.z_high
        if clamp is not None:
            ring = math.sqrt(self.radius ** 2 - clamp ** 2)
            planar = math.hypot(d[0], d[1]) or 1.0
            u = [d[0] / planar * ring, d[1] / planar * ring, clamp]
        return [self.centre[i] + u[i] for i in range(3)]


class PlaneFace:
    """Any planar face; its closest point to anything is a point inside the gap."""

    def __init__(self, anchor=(0.0, 0.0, 0.0025)):
        self.anchor = anchor

    def closest_point_trimmed(self, point):
        return list(self.anchor)


class Box:
    def __init__(self, box):
        self.box = box

    def bounding_box(self):
        return list(self.box)


class FakeCubit:
    def __init__(self, volumes, faces):
        self.volumes = volumes          # id -> (bbox, [face ids])
        self.faces = faces              # id -> face object

    def get_relatives(self, kind, entity, target):
        assert (kind, target) == ("volume", "surface")
        return list(self.volumes[entity][1])

    def get_surface_type(self, face):
        return "sphere surface" if isinstance(self.faces[face], SphereFace) else "plane surface"

    def surface(self, face):
        return self.faces[face]

    def volume(self, volume):
        return Box(self.volumes[volume][0])


def _layered_air():
    """Upper and lower air remainders plus two gap slabs, as on the C-type mesh."""
    faces = {
        1: SphereFace(ORIGIN, z_low=0.0), 2: PlaneFace(), 3: PlaneFace(),
        4: SphereFace(ORIGIN, z_high=0.0), 5: PlaneFace(), 6: PlaneFace(),
        7: PlaneFace(), 8: PlaneFace(), 9: PlaneFace(), 10: PlaneFace(),
        # an unrelated sphere (e.g. a ball inside the air) of another radius
        11: SphereFace(ORIGIN, radius=0.01),
    }
    volumes = {
        100: ((-R, -R, 0.0, R, R, R), [1, 2, 3, 11]),
        101: ((-R, -R, -R, R, R, 0.0), [4, 5, 6]),
        102: ((-0.025, -0.02, 0.0, 0.025, 0.02, 0.0025), [7, 8]),
        103: ((-0.025, -0.02, -0.0025, 0.025, 0.02, 0.0), [9, 10]),
    }
    return FakeCubit(volumes, faces)


def test_gap_slab_and_foreign_sphere_faces_are_not_kelvin_int():
    cubit = _layered_air()
    faces = add_kelvin._air_sphere_faces(cubit, [100, 101, 102, 103], R, "air")
    assert sorted(faces) == [1, 4]


def test_air_block_without_a_face_on_the_sphere_is_rejected():
    cubit = FakeCubit({200: ((-0.1, -0.1, -0.1, 0.1, 0.1, 0.1), [1, 2])},
                      {1: PlaneFace(), 2: PlaneFace()})
    with pytest.raises(RuntimeError, match="No face"):
        add_kelvin._air_sphere_faces(cubit, [200], R, "air")


def test_air_sphere_off_the_origin_is_rejected_not_half_supported():
    shifted = (0.01, 0.0, 0.0)
    cubit = FakeCubit({300: ((-R + 0.01, -R, -R, R + 0.01, R, R), [1])},
                      {1: SphereFace(shifted)})
    with pytest.raises(NotImplementedError, match="origin"):
        add_kelvin._air_sphere_faces(cubit, [300], R, "air")


def _pairing_cubit(kelvin_faces):
    faces = {1: SphereFace(ORIGIN, z_low=0.0), 4: SphereFace(ORIGIN, z_high=0.0)}
    faces.update(kelvin_faces)
    return FakeCubit({}, faces)


def test_reflected_pairing_copies_the_upper_face_and_reflects_the_lower():
    cubit = _pairing_cubit({50: SphereFace(OFFSET, z_low=0.0)})
    pairs, reflected = add_kelvin._pair_sphere_faces(
        cubit, [1, 4], [50], OFFSET, R, reflected_z=True)
    assert pairs == [(1, 50)] and reflected == [4]


def test_reflected_pairing_rejects_a_kelvin_face_that_does_not_match():
    cubit = _pairing_cubit({50: SphereFace(OFFSET, z_high=0.0)})
    with pytest.raises(RuntimeError):
        add_kelvin._pair_sphere_faces(cubit, [1, 4], [50], OFFSET, R, reflected_z=True)


def test_reflected_pairing_rejects_an_exterior_cut_off_its_centre():
    # the exterior sphere cut on a plane that does not pass through its centre
    cubit = _pairing_cubit({50: SphereFace(OFFSET, z_low=0.05)})
    with pytest.raises(RuntimeError):
        add_kelvin._pair_sphere_faces(cubit, [1, 4], [50], OFFSET, R, reflected_z=True)


def test_reflected_pairing_rejects_a_missing_lower_face():
    cubit = _pairing_cubit({50: SphereFace(OFFSET, z_low=0.0)})
    with pytest.raises(RuntimeError):
        add_kelvin._pair_sphere_faces(cubit, [1], [50], OFFSET, R, reflected_z=True)


def test_unreflected_pairing_needs_a_partner_for_every_face():
    cubit = FakeCubit({}, {1: SphereFace(ORIGIN), 50: SphereFace(OFFSET)})
    pairs, reflected = add_kelvin._pair_sphere_faces(
        cubit, [1], [50], OFFSET, R, reflected_z=False)
    assert pairs == [(1, 50)] and reflected == []
    with pytest.raises(RuntimeError):
        add_kelvin._pair_sphere_faces(cubit, [1], [], OFFSET, R, reflected_z=False)
