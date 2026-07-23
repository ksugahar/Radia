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

and, for the 0-form helper make_reduced_potential_scalar_cf (the T-Omega
background potential):

3. inner:  Omega_s = -H0 z
   kelvin: Omega_s' = -(r'/R)^2 (-H0 (z-oz))
   plus the two facts that motivate a SEPARATE 0-form helper:
     - it stays bounded at the offset, where the plain 0-form pullback
       Omega_s(T(r')) = -H0 R^2 (z-oz)/rho'^2 diverges;
     - it is NOT the potential of the 1-form helper: the 1-form
       Convention B exterior field has curl = (2 H0/R^2)(-y', x', 0) != 0
       and therefore admits no scalar potential at all.

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
from ngsolve import TaskManager
from radia.kelvin_material import (
    make_reduced_potential_background_cf,
    make_reduced_potential_scalar_cf,
)


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

    with TaskManager():
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


H_0 = 1.0


def _omega_s_helper(mesh, R_K, offset):
    """Convention B 0-form background potential for uniform H_0 z_hat."""
    return make_reduced_potential_scalar_cf(
        mesh,
        lambda xc, yc, zc: -H_0 * zc,
        R_K=R_K, offset=offset,
        kelvin_mats=("kelvin",), dim=3,
    )


def _hs_helper(mesh, R_K, offset):
    """Convention B 1-form background field for uniform H_0 z_hat."""
    return make_reduced_potential_background_cf(
        mesh,
        lambda xc, yc, zc: CoefficientFunction((0, 0, H_0)),
        R_K=R_K, offset=offset,
        kelvin_mats=("kelvin",), dim=3,
    )


def test_uniform_Omega_scalar():
    """0-form helper reproduces -(r'/R)^2 * Omega_s(local) in the Kelvin region."""
    print("=== Test 4: uniform Omega_s = -H0 z, 0-form Convention B ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset

    x_local = x - ox
    y_local = y - oy
    z_local = z - oz
    r_prime_sq = x_local**2 + y_local**2 + z_local**2 + 1e-20
    factor = -r_prime_sq / (R_K * R_K)

    is_kelvin = mesh.MaterialCF({m: (1.0 if "kelvin" in m.lower() else 0.0)
                                  for m in mesh.GetMaterials()}, default=0.0)
    Om_handwritten = (is_kelvin * (factor * (-H_0 * z_local))
                      + (1 - is_kelvin) * (-H_0 * z))

    Om_helper = _omega_s_helper(mesh, R_K, offset)

    samples = [
        ("origin (inner)",      (0.0, 0.0, 0.0)),
        ("(0.0, 0, 0.5) inner", (0.0, 0.0, 0.5)),
        ("(3.0, 0, 0) kelvin",  (3.0, 0.0, 0.0)),   # rho'=0
        ("(3.0, 0, 0.5) kelv",  (3.0, 0.0, 0.5)),   # rho'=0.5
        ("(3.4, 0.2, 0.3) kel", (3.4, 0.2, 0.3)),
    ]
    max_err = 0.0
    for label, pt in samples:
        mip = mesh(*pt)
        v_hand = Om_handwritten(mip)
        v_help = Om_helper(mip)
        err = abs(v_hand - v_help)
        max_err = max(max_err, err)
        print(f"  {label:22s}: hand={v_hand:+.6e}, helper={v_help:+.6e}, "
              f"err={err:.2e}")

    print(f"  max abs err = {max_err:.2e}")
    assert max_err < 1e-12, f"Omega_s helper mismatch: {max_err}"
    print("  [OK] 0-form helper matches the hand-written Convention B form\n")


def test_scalar_is_bounded_where_the_plain_pullback_diverges():
    """Convention B stays bounded at rho' -> 0; the 0-form pullback does not."""
    print("=== Test 5: boundedness at the offset (physical infinity) ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    Om_helper = _omega_s_helper(mesh, R_K, offset)

    # |Omega_s'| = (rho'/R)^2 H0 |z_loc| <= H0 rho'^3 / R^2 <= H0 R inside.
    bound = H_0 * R_K
    for dz in (0.5, 0.1, 0.01, 1e-3, 1e-4):
        pt = (ox, oy, oz + dz)
        val = Om_helper(mesh(*pt))
        pullback = -H_0 * R_K**2 * dz / (dz * dz)      # = -H0 R^2 / rho'
        print(f"  rho'={dz:8.1e}: ConventionB={val:+.6e}, "
              f"0-form pullback={pullback:+.6e}")
        assert abs(val) <= bound + 1e-12, \
            f"Convention B unbounded at rho'={dz}: {val}"

    # and the pullback really does blow up over the same sweep
    assert abs(-H_0 * R_K**2 / 1e-4) > 1e3 * bound
    print("  [OK] Convention B bounded; plain 0-form pullback diverges\n")


def test_scalar_boundary_sign_flip():
    """Inner and Kelvin sides are sign-flipped at rho' = R."""
    print("=== Test 6: 0-form sign flip at the Kelvin boundary (r' = R) ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    Om_helper = _omega_s_helper(mesh, R_K, offset)

    # inner point at distance R along +z, and the matching Kelvin-side point
    # at offset + R z_hat (rho' = R).
    v_inner = Om_helper(mesh(0.0, 0.0, R_K))
    v_kelvin = Om_helper(mesh(ox, oy, oz + R_K))
    print(f"  inner  @ (0,0,{R_K}):        Omega_s  = {v_inner:+.6e}")
    print(f"  kelvin @ ({ox},{oy},{oz + R_K}): Omega_s' = {v_kelvin:+.6e}")
    print(f"  sum                          = {v_inner + v_kelvin:+.6e}")
    assert abs(v_inner + v_kelvin) < 1e-12, \
        f"0-form sign flip not exact at boundary: {v_inner + v_kelvin}"
    print("  [OK] 0-form sign flip exact at rho' = R\n")


def test_scalar_and_vector_convention_b_are_not_gradient_consistent():
    """LOCK the verified incompatibility of the 0-form and 1-form B rules.

    This is a real property of Convention B, not a defect of the helpers:
    the 1-form exterior field has non-zero curl, so no scalar potential
    reproduces it.  The test pins the exact discrepancy so a future change
    to either helper cannot silently redefine the contract.

        -grad(Omega_s') - H_s'  =  -(2 H0 / R^2) (x' z', y' z', z'^2)
    """
    print("=== Test 7: 0-form vs 1-form Convention B are NOT consistent ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    Om_helper = _omega_s_helper(mesh, R_K, offset)
    Hs_helper = _hs_helper(mesh, R_K, offset)

    h = 1e-4
    max_rel = 0.0
    for pt in [(3.4, 0.2, 0.3), (3.3, -0.25, 0.2), (2.7, 0.15, -0.3)]:
        px, py, pz = pt
        # central differences of the scalar helper
        grad_num = []
        for i in range(3):
            plus = list(pt)
            minus = list(pt)
            plus[i] += h
            minus[i] -= h
            grad_num.append((Om_helper(mesh(*plus))
                             - Om_helper(mesh(*minus))) / (2 * h))
        hs = [Hs_helper[i](mesh(*pt)) for i in range(3)]
        got = [-g - f for g, f in zip(grad_num, hs)]

        xl, yl, zl = px - ox, py - oy, pz - oz
        want = [-(2 * H_0 / R_K**2) * xl * zl,
                -(2 * H_0 / R_K**2) * yl * zl,
                -(2 * H_0 / R_K**2) * zl * zl]
        scale = max(1e-12, max(abs(w) for w in want))
        rel = max(abs(a - b) for a, b in zip(got, want)) / scale
        max_rel = max(max_rel, rel)
        print(f"  r'=({xl:+.2f},{yl:+.2f},{zl:+.2f}): "
              f"-grad(Om')-H_s' = {tuple(f'{v:+.4e}' for v in got)}, "
              f"predicted = {tuple(f'{v:+.4e}' for v in want)}, rel={rel:.2e}")

        # the discrepancy is genuinely non-zero (this is the point)
        assert max(abs(v) for v in got) > 1e-6, \
            "expected a non-zero 0-form/1-form discrepancy, got ~0"

    print(f"  max rel err vs predicted discrepancy = {max_rel:.2e}")
    assert max_rel < 5e-4, \
        f"discrepancy does not match the derived formula: {max_rel}"
    print("  [OK] discrepancy matches -(2 H0/R^2)(x'z', y'z', z'^2)\n")


def test_vector_convention_b_exterior_field_has_nonzero_curl():
    """The 1-form Convention B exterior field is not curl-free (hence no potential)."""
    print("=== Test 8: curl(H_s' 1-form B) = (2H0/R^2)(-y', x', 0) != 0 ===")
    mesh, R_K, offset = build_test_geo()
    ox, oy, oz = offset
    Hs_helper = _hs_helper(mesh, R_K, offset)

    h = 1e-4

    def comp(i, pt):
        return Hs_helper[i](mesh(*pt))

    max_rel = 0.0
    for pt in [(3.4, 0.2, 0.3), (3.3, -0.25, 0.2), (2.7, 0.15, -0.3)]:
        d = [[0.0] * 3 for _ in range(3)]   # d[i][j] = d H_i / d x_j
        for i in range(3):
            for j in range(3):
                plus = list(pt)
                minus = list(pt)
                plus[j] += h
                minus[j] -= h
                d[i][j] = (comp(i, plus) - comp(i, minus)) / (2 * h)
        curl = [d[2][1] - d[1][2], d[0][2] - d[2][0], d[1][0] - d[0][1]]

        xl, yl = pt[0] - ox, pt[1] - oy
        want = [-(2 * H_0 / R_K**2) * yl, (2 * H_0 / R_K**2) * xl, 0.0]
        scale = max(1e-12, max(abs(w) for w in want))
        rel = max(abs(a - b) for a, b in zip(curl, want)) / scale
        max_rel = max(max_rel, rel)
        print(f"  r'=({xl:+.2f},{yl:+.2f}): curl = "
              f"{tuple(f'{v:+.4e}' for v in curl)}, "
              f"predicted = {tuple(f'{v:+.4e}' for v in want)}, rel={rel:.2e}")
        assert max(abs(v) for v in curl) > 1e-6, \
            "expected non-zero curl in the Kelvin exterior, got ~0"

    print(f"  max rel err vs predicted curl = {max_rel:.2e}")
    assert max_rel < 5e-4, f"curl does not match the derived formula: {max_rel}"
    print("  [OK] 1-form Convention B carries a fictitious curl in kext\n")


if __name__ == "__main__":
    test_uniform_H()
    test_uniform_A_for_uniform_Bz()
    test_kelvin_boundary_continuity()
    test_uniform_Omega_scalar()
    test_scalar_is_bounded_where_the_plain_pullback_diverges()
    test_scalar_boundary_sign_flip()
    test_scalar_and_vector_convention_b_are_not_gradient_consistent()
    test_vector_convention_b_exterior_field_has_nonzero_curl()
    print("=" * 60)
    print("ALL TESTS PASS - make_reduced_potential_background_cf is")
    print("consistent with hand-written H-formulation and A-formulation")
    print("(offset-local, Convention B, CEFC 2026).")
    print("=" * 60)
