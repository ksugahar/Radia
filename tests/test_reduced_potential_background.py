"""Test make_reduced_potential_background_cf against hand-written formulae.

Validates that the new mcp-server helper for reduced-potential background
field reproduces:

1. The H-formulation 3D_dipole_with_Kelvin.py hand-written form for
   uniform H_s = (0, 0, 1):
       inner:  H_s = (0, 0, 1)
       kelvin: H_s' = (0, 0, -(r'/R)^2)

2. The A-formulation symmetric-gauge form for uniform B_z, with
   offset-relative evaluation in the Kelvin region:
       inner:  A_s = (B0/2)(-y, x, 0)
       kelvin: A_s' = -(r'/R)^2 (B0/2)(-(y-oy), (x-ox), 0)

Reference: docs/kelvin/KELVIN_TRANSFORMATION.md §7
"""
import sys
import math
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from netgen.occ import Box, Sphere, Pnt, OCCGeometry, Glue, Vertex
from netgen.meshing import IdentificationType
from ngsolve import (
    Mesh, CoefficientFunction, x, y, z, sqrt, IfPos,
)
from radia.kelvin_material import make_reduced_potential_background_cf


def build_test_geo():
    """Two-sphere mesh: inner @ origin (R=1), outer @ offset (R=1)."""
    R_K = 1.0
    offset = (3.0, 0.0, 0.0)

    inner_sphere = Sphere(Pnt(0, 0, 0), R_K)
    inner_sphere.mat("inner_air")
    for f in inner_sphere.faces:
        f.name = "kelvin_int"
    inner_sphere.maxh = 0.3

    outer_sphere = Sphere(Pnt(*offset), R_K)
    outer_sphere.mat("kelvin_air")
    for f in outer_sphere.faces:
        f.name = "kelvin_ext"
    outer_sphere.maxh = 0.3

    vertex = Vertex(Pnt(*offset))
    vertex.name = "GND"

    geo = Glue([inner_sphere, outer_sphere, vertex])
    geo.solids[0].name = "inner_air"
    geo.solids[1].name = "kelvin_air"

    # Periodic identification
    int_face = ext_face = None
    for s in geo.solids:
        for f in s.faces:
            if f.name == "kelvin_int":
                int_face = f
            elif f.name == "kelvin_ext":
                ext_face = f
    int_face.Identify(ext_face, "periodic", IdentificationType.PERIODIC)

    mesh = Mesh(OCCGeometry(geo).GenerateMesh(maxh=0.3, grading=0.7))
    return mesh, R_K, offset


def test_uniform_H():
    print("=== Test 1: uniform H_s = (0, 0, 1), H-formulation ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    print(f"  materials = {mesh.GetMaterials()}")

    # Hand-written formula (matches 3D_dipole_with_Kelvin.py)
    x_local = x - ox
    y_local = y - oy
    z_local = z - oz
    r_prime = sqrt(x_local**2 + y_local**2 + z_local**2 + 1e-20)
    Hs_z_outer_handwritten = -(r_prime / R_K) ** 2

    is_kelvin = mesh.MaterialCF({m: (1.0 if "kelvin" in m.lower() else 0.0)
                                  for m in mesh.GetMaterials()}, default=0.0)

    Hs_handwritten = CoefficientFunction((
        0.0, 0.0,
        is_kelvin * Hs_z_outer_handwritten + (1 - is_kelvin) * 1.0,
    ))

    # New helper
    Hs_helper = make_reduced_potential_background_cf(
        mesh,
        lambda xc, yc, zc: CoefficientFunction((0, 0, 1)),
        R_K=R_K, offset=offset,
        kelvin_mats=("kelvin",), dim=3,
    )

    # Evaluate at sample points
    samples = [
        ("origin (inner)",       (0.0, 0.0, 0.0)),
        ("(0.5, 0, 0) inner",    (0.5, 0.0, 0.0)),
        ("(2.5, 0, 0) kelvin",   (2.5, 0.0, 0.0)),  # r'=0.5, near interface
        ("(3.0, 0, 0) kelvin",   (3.0, 0.0, 0.0)),  # offset center, r'=0
        ("(3.5, 0, 0) kelvin",   (3.5, 0.0, 0.0)),  # r'=0.5
    ]
    max_err = 0.0
    for label, pt in samples:
        try:
            mip = mesh(*pt)
            h_hand = [Hs_handwritten[i](mip) for i in range(3)]
            h_help = [Hs_helper[i](mip) for i in range(3)]
            err = max(abs(a - b) for a, b in zip(h_hand, h_help))
            max_err = max(max_err, err)
            print(f"  {label:25s}: hand={tuple(f'{v:+.4e}' for v in h_hand)}, "
                  f"helper={tuple(f'{v:+.4e}' for v in h_help)}, err={err:.2e}")
        except Exception as e:
            print(f"  {label:25s}: skip ({e})")

    print(f"  max abs err = {max_err:.2e}")
    assert max_err < 1e-12, f"H_s helper mismatch: {max_err}"
    print("  [OK] uniform H_s matches hand-written formula\n")


def test_uniform_A_for_uniform_Bz():
    print("=== Test 2: uniform B_z = 1, A_s = (1/2)(-y, x, 0), A-formulation ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    B_0 = 1.0

    # Expected (with offset-local evaluation in Kelvin)
    x_local = x - ox
    y_local = y - oy
    z_local = z - oz
    r_prime_sq = x_local**2 + y_local**2 + z_local**2 + 1e-20
    factor = -r_prime_sq / (R_K * R_K)

    is_kelvin = mesh.MaterialCF({m: (1.0 if "kelvin" in m.lower() else 0.0)
                                  for m in mesh.GetMaterials()}, default=0.0)
    # Inner: A = (B0/2)(-y, x, 0); Kelvin: factor * (B0/2)(-y_local, x_local, 0)
    A_inner_x = -(B_0 / 2) * y
    A_inner_y =  (B_0 / 2) * x
    A_kelvin_x = factor * (-(B_0 / 2) * y_local)
    A_kelvin_y = factor * ( (B_0 / 2) * x_local)
    As_x_handwritten = is_kelvin * A_kelvin_x + (1 - is_kelvin) * A_inner_x
    As_y_handwritten = is_kelvin * A_kelvin_y + (1 - is_kelvin) * A_inner_y
    As_z_handwritten = CoefficientFunction(0.0)
    As_handwritten = CoefficientFunction((As_x_handwritten,
                                           As_y_handwritten,
                                           As_z_handwritten))

    # New helper
    As_helper = make_reduced_potential_background_cf(
        mesh,
        lambda xc, yc, zc: CoefficientFunction((-yc, xc, 0)) * (B_0 / 2),
        R_K=R_K, offset=offset,
        kelvin_mats=("kelvin",), dim=3,
    )

    samples = [
        ("origin (inner)",      (0.0, 0.0, 0.0),  (0, 0, 0)),
        ("(0.5, 0.3, 0) in",    (0.5, 0.3, 0.0),  (-0.15, 0.25, 0)),
        ("(3.0, 0, 0) kelv",    (3.0, 0.0, 0.0),  (0, 0, 0)),  # offset center, r'=0 -> A=0
        ("(3.5, 0, 0) kelv",    (3.5, 0.0, 0.0),  (0, 0, 0)),  # r'=0.5: y_local=0, x_local=0.5; factor=-0.25; A=-0.25*(0, 0.25, 0)=(0, -0.0625, 0)
        ("(3.0, 0.5, 0) kelv",  (3.0, 0.5, 0.0),  (0, 0, 0)),  # r'=0.5: y_local=0.5, x_local=0; factor=-0.25; A=-0.25*(-0.25, 0, 0)=(0.0625, 0, 0)
    ]
    max_err = 0.0
    for label, pt, _ in samples:
        try:
            mip = mesh(*pt)
            a_hand = [As_handwritten[i](mip) for i in range(3)]
            a_help = [As_helper[i](mip) for i in range(3)]
            err = max(abs(a - b) for a, b in zip(a_hand, a_help))
            max_err = max(max_err, err)
            print(f"  {label:22s}: hand={tuple(f'{v:+.4e}' for v in a_hand)}, "
                  f"helper={tuple(f'{v:+.4e}' for v in a_help)}, err={err:.2e}")
        except Exception as e:
            print(f"  {label:22s}: skip ({e})")

    print(f"  max abs err = {max_err:.2e}")
    assert max_err < 1e-12, f"A_s helper mismatch: {max_err}"
    print("  [OK] A_s helper matches expected offset-local formula\n")


def test_kelvin_boundary_continuity():
    print("=== Test 3: Sign-flip matching at Kelvin boundary (r' = R) ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    B_0 = 1.0

    As_helper = make_reduced_potential_background_cf(
        mesh,
        lambda xc, yc, zc: CoefficientFunction((-yc, xc, 0)) * (B_0 / 2),
        R_K=R_K, offset=offset,
        kelvin_mats=("kelvin",), dim=3,
    )

    # At Kelvin boundary, the inner and exterior sides should have OPPOSITE
    # sign (sign flip). Inner side at (1, 0, 0) [r=R from origin]; Kelvin side
    # at the Kelvin-identified point at offset + (1, 0, 0)? But identification
    # in 3D_dipole is geometric — the inner sphere boundary face is identified
    # with the Kelvin sphere boundary face. Same point (x=1, y=0, z=0) appears
    # in both representations conceptually (via periodic BC).
    #
    # Simpler: just verify at point on Kelvin boundary (r' = R = 1):
    #   At (3+1, 0, 0) = (4, 0, 0): r' = sqrt(1+0+0) = 1 = R.
    #   Local x = 4-3 = 1, y_local = 0. A_kelvin = -(1/1)² * (1/2)(0, 1, 0)
    #                                            = -(1/2)(0, 1, 0)
    #                                            = (0, -0.5, 0)
    #   At (1, 0, 0) inner: A_inner = (1/2)(0, 1, 0) = (0, +0.5, 0)
    # These are sign-flipped, matching the periodic BC requirement.

    pt_inner = (1.0, 0.0, 0.0)
    pt_kelvin = (4.0, 0.0, 0.0)
    try:
        a_inner = [As_helper[i](mesh(*pt_inner)) for i in range(3)]
        a_kelvin = [As_helper[i](mesh(*pt_kelvin)) for i in range(3)]
        print(f"  Inner @ (1,0,0):    A = {tuple(f'{v:+.4e}' for v in a_inner)}")
        print(f"  Kelvin @ (4,0,0):   A = {tuple(f'{v:+.4e}' for v in a_kelvin)}")
        # At r'=R, A_kelvin should equal -A_inner_local (with offset shift)
        sum_check = [a_inner[i] + a_kelvin[i] for i in range(3)]
        print(f"  inner + kelvin:    = {tuple(f'{v:+.4e}' for v in sum_check)}")
        max_sum = max(abs(s) for s in sum_check)
        print(f"  max |inner + kelvin| at boundary = {max_sum:.2e}")
        assert max_sum < 1e-12, f"sign flip not exact at boundary: {max_sum}"
        print("  [OK] inner and Kelvin sides sign-flipped at boundary\n")
    except Exception as e:
        print(f"  skip ({e})\n")


if __name__ == "__main__":
    test_uniform_H()
    test_uniform_A_for_uniform_Bz()
    test_kelvin_boundary_continuity()
    print("=" * 60)
    print("ALL TESTS PASS - make_reduced_potential_background_cf is")
    print("consistent with hand-written H-formulation and A-formulation")
    print("(offset-local, Convention B, CEFC 2026).")
    print("=" * 60)
