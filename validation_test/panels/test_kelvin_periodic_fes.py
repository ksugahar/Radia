"""FES-only verification of Kelvin Periodic identification.

Sub-second tests (no physics solve) that lock the 2026-04-25 finding:
both OCC `Identify(IdentificationType.PERIODIC)` and Cubit
`export netgen` produce .vol files where NGSolve's `Periodic`
H1 FES correctly slaves high-order DOFs at order=2.

Catches regression of either:
  * `detect_kelvin_offset` (must return mean(kelvin_ext) -
    mean(kelvin_int), NOT the kelvin region centroid).
  * The Periodic FES coupling itself (slaved DOFs and
    Set(1)|kelvin_int -> kelvin_ext boundary integral ratio).
"""
import os
import sys
import pytest

PANELS = os.path.join(os.path.dirname(__file__), "..", "..", "src", "radia", "panels")
sys.path.insert(0, PANELS)


def _build_occ_kelvin_mesh(maxh_air=0.08):
    """Build the canonical OCC Kelvin pattern in-memory:
    `Identify(PERIODIC)` AFTER `Glue([...])`.  Returns the NGSolve mesh."""
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


def test_occ_kelvin_p1_slaves_dofs():
    """OCC Identify-after-Glue: Periodic at p=1 slaves > 0 DOFs."""
    from ngsolve import H1, Periodic
    mesh = _build_occ_kelvin_mesh()
    fb = H1(mesh, order=1, dirichlet="GND")
    fp = Periodic(fb)
    slaved = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
    assert slaved > 0, f"p=1 Periodic slaved 0 DOFs (got slaved={slaved})"


def test_occ_kelvin_p2_slaves_more_than_p1():
    """OCC Identify-after-Glue: Periodic at p=2 slaves ~4x p=1.

    Key claim: high-order edge / face DOFs ARE coupled by Periodic
    FES even when NGSolve's identifications API only exposes point
    pairs.  If high-order coupling broke (regression), p=2 slaved
    would equal p=1 slaved.
    """
    from ngsolve import H1, Periodic
    mesh = _build_occ_kelvin_mesh()
    with TaskManager():
        mesh.Curve(2)
        slaved_at = {}
        for p in (1, 2):
            fb = H1(mesh, order=p, dirichlet="GND")
            fp = Periodic(fb)
            slaved_at[p] = sum(fb.FreeDofs()) - sum(fp.FreeDofs())
        assert slaved_at[2] > slaved_at[1] * 2, (
            f"p=2 should slave > 2x p=1 DOFs (got p1={slaved_at[1]}, "
            f"p2={slaved_at[2]})")


def test_occ_kelvin_set_then_integrate_ratio():
    """OCC Identify-after-Glue: Set(1) on kelvin_int MUST equal kelvin_ext."""
    from ngsolve import H1, Periodic, GridFunction, Integrate
    from ngsolve import TaskManager
    mesh = _build_occ_kelvin_mesh()
    mesh.Curve(2)
    fp = Periodic(H1(mesh, order=2, dirichlet="GND"))
    gfu = GridFunction(fp)
    gfu.vec[:] = 0
    gfu.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
    a_int = float(Integrate(gfu * gfu, mesh,
                            definedon=mesh.Boundaries("kelvin_int")))
    a_ext = float(Integrate(gfu * gfu, mesh,
                            definedon=mesh.Boundaries("kelvin_ext")))
    assert a_int > 0, "kelvin_int integral is zero"
    ratio = a_ext / a_int
    assert abs(ratio - 1.0) < 1e-3, (
        f"Periodic Kelvin ratio not 1.0: int={a_int:.4e}, "
        f"ext={a_ext:.4e}, ratio={ratio:.4f}")


def test_detect_kelvin_offset_returns_translation_not_centroid():
    """detect_kelvin_offset must return mean(kelvin_ext) - mean(kelvin_int).

    Regression guard: pre-2026-04-25 implementation returned the kelvin
    REGION centroid (i.e., the average of all kelvin-material vertex
    positions), which for a partial-octant model gave a translation
    that was wrong by ~0.146 m in x and z.
    """
    if "calc_common" in sys.modules:
        del sys.modules["calc_common"]
    from calc_common import detect_kelvin_offset

    mesh = _build_occ_kelvin_mesh()
    offset = detect_kelvin_offset(mesh)
    # Expected exact translation: (0, 0, 3.0)
    assert abs(offset[0]) < 1e-6, f"x-component nonzero: {offset[0]}"
    assert abs(offset[1]) < 1e-6, f"y-component nonzero: {offset[1]}"
    assert abs(offset[2] - 3.0) < 1e-3, f"z-component != 3.0: {offset[2]}"


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v"]))
