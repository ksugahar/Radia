"""Test scalar BIE + SIBC solver against analytical sphere solution."""

import math
import sys
import os
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.join(ROOT, "src", "radia"))

MU_0 = 4e-7 * np.pi


def test_sphere_pec():
    """Scalar BIE on PEC sphere should match analytical (3/2)*H0*sqrt(2/3)."""
    from ngsolve import Mesh, CF, z
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    from scalar_bie_sibc import ScalarBIE_SIBC

    R = 0.01
    B0 = 0.001
    H0 = B0 / MU_0
    freq = 1000
    omega = 2 * math.pi * freq

    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=R / 3))
    with TaskManager():
        mesh.Curve(2)

        solver = ScalarBIE_SIBC(mesh, order=1)
        print(f"Assembly: {solver.t_assembly:.1f}s, {solver.ndof} DOFs")

        # PEC: Z_s = 0
        H_inc = CF((0, 0, H0))
        phi_inc = CF(-H0 * z)
        result = solver.solve(H_inc, Zs=0, omega=omega, phi_inc_cf=phi_inc)

        H_ana = 1.5 * H0 * math.sqrt(2.0 / 3.0)
        err = abs(result['H_rms'] / H_ana - 1)
        print(f"PEC: H_rms = {result['H_rms']:.2f}, analytical = {H_ana:.2f}, "
              f"error = {err*100:.2f}%")
        assert err < 0.01, f"PEC error {err*100:.2f}% > 1%"


def test_sphere_sibc_sweep():
    """Scalar BIE + SIBC should match analytical for all Z_s ratios."""
    from ngsolve import Mesh, CF, z
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    from scalar_bie_sibc import ScalarBIE_SIBC

    R = 0.01
    B0 = 0.001
    H0 = B0 / MU_0
    freq = 1000
    omega = 2 * math.pi * freq
    beta = 1j * omega * MU_0 * R

    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=R / 3))
    mesh.Curve(2)

    solver = ScalarBIE_SIBC(mesh, order=1)

    H_inc = CF((0, 0, H0))
    phi_inc = CF(-H0 * z)

    print(f"\n{'ratio':>8s} {'Ana':>10s} {'BIE':>10s} {'err%':>8s}")
    print("-" * 40)

    max_err = 0
    for ratio in [0.001, 0.1, 1.0, 5.0, 10.0]:
        Zs = ratio * abs(beta) * (1 + 1j) / math.sqrt(2)
        J_max_ana = abs(1.5 * 1j * omega * B0 * R / (beta + Zs))
        H_ana = J_max_ana * math.sqrt(2.0 / 3.0)

        result = solver.solve(H_inc, Zs, omega, phi_inc_cf=phi_inc)
        err = abs(result['H_rms'] / H_ana - 1)
        max_err = max(max_err, err)
        print(f"{ratio:8.3f} {H_ana:10.2f} {result['H_rms']:10.2f} "
              f"{err*100:+8.2f}%")

    print(f"\nMax error: {max_err*100:.2f}%")
    assert max_err < 0.01, f"Max error {max_err*100:.2f}% > 1%"


def test_phi_inc_from_H():
    """Test phi_inc reconstruction from H_inc via surface Poisson."""
    from ngsolve import Mesh, CF, z
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    from scalar_bie_sibc import ScalarBIE_SIBC

    R = 0.01
    B0 = 0.001
    H0 = B0 / MU_0
    freq = 1000
    omega = 2 * math.pi * freq
    beta = 1j * omega * MU_0 * R

    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=R / 3))
    mesh.Curve(2)

    solver = ScalarBIE_SIBC(mesh, order=1)

    H_inc = CF((0, 0, H0))
    phi_inc_direct = CF(-H0 * z)

    # Compare: direct phi_inc vs reconstructed from H_inc
    Zs = 1.0 * abs(beta) * (1 + 1j) / math.sqrt(2)

    r_direct = solver.solve(H_inc, Zs, omega, phi_inc_cf=phi_inc_direct)
    r_recon = solver.solve(H_inc, Zs, omega, phi_inc_cf=None)

    err = abs(r_direct['H_rms'] / r_recon['H_rms'] - 1)
    print(f"\nDirect phi_inc: H_rms = {r_direct['H_rms']:.4f}")
    print(f"Reconstructed:  H_rms = {r_recon['H_rms']:.4f}")
    print(f"Difference: {err*100:.2f}%")
    assert err < 0.02, f"phi_inc reconstruction error {err*100:.2f}% > 2%"


def test_frequency_sweep():
    """Test frequency sweep API."""
    from ngsolve import Mesh, CF, z
    from ngsolve import TaskManager
    from netgen.occ import Sphere, Pnt, OCCGeometry, Glue
    from scalar_bie_sibc import ScalarBIE_SIBC

    R = 0.01
    B0 = 0.001
    H0 = B0 / MU_0
    sigma_cu = 5.8e7  # copper

    sph = Sphere(Pnt(0, 0, 0), R)
    for f in sph.faces:
        f.name = "sphere"
    geo = OCCGeometry(Glue(sph.faces))
    mesh = Mesh(geo.GenerateMesh(maxh=R / 3))
    mesh.Curve(2)

    solver = ScalarBIE_SIBC(mesh, order=1)
    phi_inc = CF(-H0 * z)
    H_inc = CF((0, 0, H0))

    results = solver.frequency_sweep(H_inc, freqs=[1e3, 10e3, 100e3],
                                      sigma=sigma_cu, phi_inc_cf=phi_inc)

    print(f"\n{'freq':>10s} {'delta_mm':>10s} {'H_rms':>10s} {'P_loss':>12s}")
    print("-" * 46)
    for r in results:
        print(f"{r['frequency']:10.0f} {r['skin_depth']*1e3:10.3f} "
              f"{r['H_rms']:10.2f} {r['P_loss']:12.6e}")

    # H_rms should INCREASE with frequency for copper sphere:
    # higher f -> smaller skin depth -> Z_s -> 0 (PEC limit) -> J -> J_PEC
    assert results[0]['H_rms'] < results[-1]['H_rms'], \
        "H_rms should increase with frequency (approaching PEC limit)"
    # P_loss should also increase (more current, higher R_s)
    assert results[0]['P_loss'] < results[-1]['P_loss'], \
        "P_loss should increase with frequency"
    print("\nFrequency sweep: H_rms increases toward PEC limit (OK)")


if __name__ == "__main__":
    test_sphere_pec()
    test_sphere_sibc_sweep()
    test_phi_inc_from_H()
    test_frequency_sweep()
    print("\n=== All tests passed ===")
