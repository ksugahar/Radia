"""OCC-produced periodic structure and existing H1 validation.

Assumes: producer-side OCC Identify before meshing, with unchanged geometry and
existing numerical assertions. Consumer verification never repairs the mesh.
"""
from __future__ import annotations

import pytest


def _build_occ_kelvin_mesh(maxh_air: float = 0.08):
    """Build the canonical 2-sphere Kelvin geometry with OCC Identify.

    Identify is called BEFORE GenerateMesh so the resulting mesh has
    1:1 kelvin_int / kelvin_ext vertex correspondence (the only way
    to get this in pure NGSolve without Cubit).  Returns the NGSolve
    Mesh with identifications already populated.
    """
    from netgen.occ import (Sphere, Pnt, OCCGeometry, IdentificationType,
                            Vertex, Glue)
    from ngsolve import Mesh

    sphere_radius = 0.5
    kelvin_radius = 1.0
    offset_z = 3.0

    mag_sphere = Sphere(Pnt(0, 0, 0), sphere_radius)
    mag_sphere.mat("magnetic")
    mag_sphere.maxh = 0.05
    for face in mag_sphere.faces:
        face.name = "sphere"

    inner_sphere = Sphere(Pnt(0, 0, 0), kelvin_radius)
    inner_sphere.maxh = maxh_air
    for face in inner_sphere.faces:
        face.name = "kelvin_int"
    inner_air = inner_sphere - mag_sphere
    inner_air.mat("air_inner")

    outer_sphere = Sphere(Pnt(0, 0, offset_z), kelvin_radius)
    outer_sphere.maxh = maxh_air
    outer_sphere.mat("air_outer")
    for face in outer_sphere.faces:
        face.name = "kelvin_ext"

    vertex = Vertex(Pnt(0, 0, offset_z))
    vertex.name = "GND"
    geo = Glue([inner_air, mag_sphere, outer_sphere, vertex])
    geo.solids[0].name = "air_inner"
    geo.solids[1].name = "magnetic"
    geo.solids[2].name = "air_outer"

    kint = next(f for s in geo.solids for f in s.faces if f.name == "kelvin_int")
    kext = next(f for s in geo.solids for f in s.faces if f.name == "kelvin_ext")
    kint.Identify(kext, "periodic", IdentificationType.PERIODIC)

    return Mesh(OCCGeometry(geo).GenerateMesh(maxh=maxh_air, grading=0.5))


# ----------------------------------------------------------------------
# detect_kelvin_offset
# ----------------------------------------------------------------------

def test_detect_kelvin_offset_returns_translation():
    """Auto-detected offset = mean(kelvin_ext) - mean(kelvin_int)."""
    from radia.kelvin_identify_ngsolve import detect_kelvin_offset
    mesh = _build_occ_kelvin_mesh()
    offset = detect_kelvin_offset(mesh)
    # Synthetic offset = (0, 0, 3); sub-mm jitter from Netgen surface mesh
    # vertex distribution is tolerated.
    assert abs(offset[0]) < 1e-3, f"x != 0: {offset[0]}"
    assert abs(offset[1]) < 1e-3, f"y != 0: {offset[1]}"
    assert abs(offset[2] - 3.0) < 1e-3, f"z != 3.0: {offset[2]}"


def test_detect_kelvin_offset_returns_zero_when_labels_missing():
    """Missing kelvin_int / kelvin_ext -> (0, 0, 0) (NOT an exception)."""
    from radia.kelvin_identify_ngsolve import detect_kelvin_offset
    from ngsolve import Mesh
    from netgen.occ import Sphere, Pnt, OCCGeometry
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = Mesh(geo.GenerateMesh(maxh=0.3))
    assert detect_kelvin_offset(mesh) == (0.0, 0.0, 0.0)


# ----------------------------------------------------------------------
# read-only structural verification
# ----------------------------------------------------------------------







# ----------------------------------------------------------------------
# retired repair API removed
# ----------------------------------------------------------------------











# ----------------------------------------------------------------------
# Periodic FES coupling at order 1 and 2
# ----------------------------------------------------------------------

def test_periodic_fes_couples_at_order_1():
    """Periodic FES on the OCC-Identify mesh slaves > 0 DOFs at p=1.

    This is the existing positive control from test_kelvin_periodic_fes;
    we replicate it here as a smoke check that the helper module
    integration doesn't break the underlying mesh.
    """
    from ngsolve import H1, Periodic
    mesh = _build_occ_kelvin_mesh()
    fb = H1(mesh, order=1, dirichlet="GND")
    fp = Periodic(fb)
    slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
    assert slaved > 0, f"p=1 Periodic slaved 0 DOFs (got {slaved})"


def test_periodic_set_then_integrate_ratio():
    """Set(1) on kelvin_int -> Integral(kelvin_ext) ratio = 1.0."""
    from ngsolve import H1, Periodic, GridFunction, Integrate
    from ngsolve import TaskManager
    mesh = _build_occ_kelvin_mesh()
    with TaskManager():
        mesh.Curve(2)
        fp = Periodic(H1(mesh, order=2, dirichlet="GND"))
        gfu = GridFunction(fp)
        gfu.vec[:] = 0
        gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
        a_int = float(Integrate(gfu * gfu, mesh,
                                definedon=mesh.Boundaries("kelvin_int")))
        a_ext = float(Integrate(gfu * gfu, mesh,
                                definedon=mesh.Boundaries("kelvin_ext")))
        assert a_int > 0
        ratio = a_ext / a_int
        assert abs(ratio - 1.0) < 1e-3, (
            f"Periodic Kelvin ratio not 1.0: int={a_int:.4e}, "
            f"ext={a_ext:.4e}, ratio={ratio:.4f}")


if __name__ == "__main__":
    import sys
    sys.exit(pytest.main([__file__, "-v"]))


def test_verify_occ_producer_is_read_only():
    from radia.kelvin_identify_ngsolve import verify_kelvin_identification
    mesh = _build_occ_kelvin_mesh()
    before = tuple(mesh.ngmesh.GetIdentifications())
    info = verify_kelvin_identification(mesh, master="kelvin_int", slave="kelvin_ext", group=1)
    assert info["verified"] is True
    assert info["n_pairs"] == len(before)
    assert info["name_status"] == "NOT_VERIFIABLE_FROM_THIS_API"
    assert tuple(mesh.ngmesh.GetIdentifications()) == before
