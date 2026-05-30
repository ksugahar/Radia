"""Post-hoc Kelvin identification on an existing mesh (no Cubit / OCC Identify).

Tests the public API in `radia.kelvin_identify_ngsolve`.

Pre-condition note
------------------
The post-hoc helper REQUIRES that kelvin_int and kelvin_ext surfaces
have a 1:1 vertex correspondence (modulo the Kelvin rigid translation).
Two upstream paths produce such meshes:
  - Cubit `copy mesh source vertex ... target vertex` (used by
    `add_kelvin_cubit`)
  - NGSolve OCC `Identify(IdentificationType.PERIODIC)` called BEFORE
    `GenerateMesh()`

A mesh where the two surfaces are meshed INDEPENDENTLY (e.g. two
separately defined spheres glued without Identify) is NOT supported --
Netgen produces different vertex distributions on each surface, the
Hungarian pairing gives garbage pairs, and the resulting Periodic FES
does not couple correctly.  Verified 2026-05-13 by trying it: a
no-Identify mesh gave max_dist=5.8cm with most pairs at 5cm separation
and slaved=0 DOFs in Periodic FES.

These tests therefore exercise the helper against an OCC-Identify-built
mesh (1:1 by construction).  The `skip_if_existing=False` flag forces
the helper to re-pair vertices that already have identifications from
OCC -- the same code path that Cubit-1/8-octant rescue uses in practice.
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
# has_kelvin_identification
# ----------------------------------------------------------------------

def test_has_kelvin_identification_on_occ_identify_mesh():
    """Mesh built WITH OCC Identify must report identifications present."""
    from radia.kelvin_identify_ngsolve import has_kelvin_identification
    mesh = _build_occ_kelvin_mesh()
    assert has_kelvin_identification(mesh) is True


def test_has_kelvin_identification_on_unidentified_mesh():
    """A plain sphere mesh has no identifications."""
    from radia.kelvin_identify_ngsolve import has_kelvin_identification
    from ngsolve import Mesh
    from netgen.occ import Sphere, Pnt, OCCGeometry
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = Mesh(geo.GenerateMesh(maxh=0.3))
    assert has_kelvin_identification(mesh) is False


# ----------------------------------------------------------------------
# add_kelvin_identification: API contract
# ----------------------------------------------------------------------

def test_add_kelvin_identification_skip_if_existing_default():
    """When ngmesh already has identifications, the helper must skip."""
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    mesh = _build_occ_kelvin_mesh()
    info = add_kelvin_identification(mesh)   # skip_if_existing=True default
    assert info["was_pre_existing"] is True
    assert info["n_pairs"] == 0


def test_add_kelvin_identification_force_re_add_pairs_vertices():
    """skip_if_existing=False forces re-pair on a 1:1 (OCC-Identify) mesh.

    The post-hoc helper requires a 1:1 mesh (see module docstring).
    OCC Identify before GenerateMesh produces such a mesh.  On this
    input the Hungarian pairing should:
      - register a non-zero number of pairs
      - have ZERO unmatched vertices (full 1:1 coverage)
      - have a tiny max_dist (sub-mm; Netgen surface mesh jitter)
    """
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    mesh = _build_occ_kelvin_mesh()
    info = add_kelvin_identification(mesh, skip_if_existing=False)
    assert info["was_pre_existing"] is False
    assert info["n_pairs"] > 0, f"Expected > 0 pairs, got {info}"
    assert info["n_unmatched"] == 0, (
        f"OCC-Identify mesh should have ZERO unmatched vertices "
        f"(1:1 mesh by construction).  Got: {info}")
    # max_dist should be tiny: the surfaces are identical up to numerical
    # representation (sub-mm rounding from Netgen vertex storage).
    assert info["max_dist"] < 1e-3, (
        f"max_dist too large on 1:1 OCC-Identify mesh: "
        f"{info['max_dist']:.3e}.  Expected sub-mm.  Info: {info}")


def test_add_kelvin_identification_raises_on_missing_labels():
    """A mesh without kelvin_int / kelvin_ext must raise ValueError, NOT
    silently no-op (fail-fast policy, CLAUDE.md No Fallbacks)."""
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    from ngsolve import Mesh
    from netgen.occ import Sphere, Pnt, OCCGeometry
    geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
    mesh = Mesh(geo.GenerateMesh(maxh=0.3))
    with pytest.raises(ValueError, match="kelvin"):
        add_kelvin_identification(mesh)


def test_add_kelvin_identification_explicit_offset_overrides_auto_detect():
    """When kelvin_offset is provided, auto-detect is bypassed."""
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    mesh = _build_occ_kelvin_mesh()
    info = add_kelvin_identification(
        mesh, kelvin_offset=(0.0, 0.0, 3.0), skip_if_existing=False)
    assert info["kelvin_offset"] == (0.0, 0.0, 3.0)
    assert info["n_pairs"] > 0


def test_add_kelvin_identification_wrong_offset_leaves_pairs_unmatched():
    """An obviously-wrong offset should produce all-unmatched vertices,
    not silently mis-pair.  Tests the tolerance gate."""
    from radia.kelvin_identify_ngsolve import add_kelvin_identification
    mesh = _build_occ_kelvin_mesh()
    info = add_kelvin_identification(
        mesh,
        kelvin_offset=(100.0, 0.0, 0.0),   # nonsense: 100 m away
        point_tolerance=0.01,                # tight tolerance
        skip_if_existing=False)
    assert info["n_pairs"] == 0, (
        f"Wrong offset should produce zero matched pairs at 1 cm "
        f"tolerance.  Got: {info}")
    assert info["n_unmatched"] > 0


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
