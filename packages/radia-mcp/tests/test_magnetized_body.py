r"""Uniformly-magnetized bodies -- closed-form regression tests.

Each test pins one analytic function to its exact textbook value AND an
independent, non-tautological gate (a limiting case that reproduces a DIFFERENT
known result, a special value, or a scaling law), plus the ValueError guard.
References: Jackson Sec. 5.10/5.11; Griffiths; Furlani Ch. 3; Aharoni,
J. Appl. Phys. 83, 3432 (1998); Joseph, J. Appl. Phys. 37, 4639 (1966);
Halbach, NIM 169, 1 (1980); Mallinson (1973).
"""
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.magnetized_body import (
    MU0,
    coaxial_magnets_axial_force,
    coaxial_magnets_force_gap_sweep_summary,
    cylinder_central_demag_factor,
    halbach_segmentation_factor,
    magnetic_polarizability_sphere,
    rectangular_prism_demag_factors,
    uniform_magnetized_annulus_axial_field,
    uniform_magnetized_cylinder_axial_field,
    uniform_magnetized_sphere_field,
)

M = 8.0e5          # magnetization [A/m]
A = 0.01           # a generic radius [m]


def test_sphere_field():
    # exact interior field: B_in = (2/3) mu0 M, H_in = -M/3 (Jackson 5.10)
    inside = uniform_magnetized_sphere_field(M, A, 0.0, 0.0, 0.0)
    assert inside["region"] == "interior"
    assert math.isclose(inside["Bz"], (2.0 / 3.0) * MU0 * M, rel_tol=1e-14)
    assert math.isclose(inside["H_in_z"], -M / 3.0, rel_tol=1e-14)
    assert math.isclose(inside["moment"], (4.0 / 3.0) * math.pi * A ** 3 * M,
                        rel_tol=1e-14)
    # GATE: exterior on-axis field is an EXACT point dipole B = mu0 m / (2 pi z^3)
    # (a different, independent closed form than the interior uniform field)
    z = 0.1
    ext = uniform_magnetized_sphere_field(M, A, 0.0, 0.0, z)
    assert ext["region"] == "exterior"
    m = (4.0 / 3.0) * math.pi * A ** 3 * M
    assert math.isclose(ext["Bz"], MU0 * m / (2.0 * math.pi * z ** 3), rel_tol=1e-12)
    assert math.isclose(ext["Bx"], 0.0, abs_tol=1e-18)
    with pytest.raises(ValueError):
        uniform_magnetized_sphere_field(M, -A, 0, 0, 0)


def test_cylinder_axial_field():
    L = 0.04
    # exact centre value Bz(0) = mu0 M L / sqrt(L^2 + 4 a^2)
    res0 = uniform_magnetized_cylinder_axial_field(M, A, L, 0.0)
    centre_exact = MU0 * M * L / math.sqrt(L * L + 4.0 * A * A)
    assert math.isclose(res0["Bz"], centre_exact, rel_tol=1e-14)
    assert math.isclose(res0["centre"], centre_exact, rel_tol=1e-14)
    assert math.isclose(res0["moment"], math.pi * A * A * L * M, rel_tol=1e-14)
    # GATE: a LONG rod (L >> a) has centre field -> mu0 M (open-circuit slug)
    long_rod = uniform_magnetized_cylinder_axial_field(M, A, 5.0, 0.0)
    assert math.isclose(long_rod["centre"], MU0 * M, rel_tol=1e-4)
    assert long_rod["centre"] < MU0 * M       # always below the ideal limit
    with pytest.raises(ValueError):
        uniform_magnetized_cylinder_axial_field(M, A, -L, 0.0)


def test_annulus_axial_field():
    L = 0.03
    ri, ro = 0.005, 0.012
    z = 0.02
    ring = uniform_magnetized_annulus_axial_field(M, ri, ro, L, z)["Bz"]
    # exact: ring = outer solid cylinder minus inner solid cylinder
    outer = uniform_magnetized_cylinder_axial_field(M, ro, L, z)["Bz"]
    inner = uniform_magnetized_cylinder_axial_field(M, ri, L, z)["Bz"]
    assert math.isclose(ring, outer - inner, rel_tol=1e-14)
    # GATE: as the bore shrinks (ri -> 0) the annulus tends to the solid cylinder
    near_solid = uniform_magnetized_annulus_axial_field(M, 1e-6, ro, L, z)["Bz"]
    assert math.isclose(near_solid, outer, rel_tol=1e-6)
    with pytest.raises(ValueError):
        uniform_magnetized_annulus_axial_field(M, ro, ri, L, z)   # ro <= ri


def test_prism_demag_factors():
    # exact cube: Nx = Ny = Nz = 1/3, sum = 1 (Aharoni 1998)
    cube = rectangular_prism_demag_factors(2.0, 2.0, 2.0)
    assert math.isclose(cube["Nx"], 1.0 / 3.0, rel_tol=1e-12)
    assert math.isclose(cube["Ny"], 1.0 / 3.0, rel_tol=1e-12)
    assert math.isclose(cube["Nz"], 1.0 / 3.0, rel_tol=1e-12)
    assert math.isclose(cube["sum"], 1.0, rel_tol=1e-12)
    # GATE: the trace sum = 1 holds for an ARBITRARY (non-cubic) prism, and a
    # plate thin along z (large a,b / small c) has Nz -> 1 (the slab limit)
    prism = rectangular_prism_demag_factors(3.0, 1.5, 0.7)
    assert math.isclose(prism["sum"], 1.0, rel_tol=1e-9)
    plate = rectangular_prism_demag_factors(50.0, 50.0, 1.0)
    assert plate["Nz"] > 0.9 and plate["Nx"] < 0.05
    with pytest.raises(ValueError):
        rectangular_prism_demag_factors(0.0, 1.0, 1.0)


def test_cylinder_central_demag():
    # exact: N = 1 - g/sqrt(g^2+1) and it ties to the cylinder centre field
    g = 1.0
    N = cylinder_central_demag_factor(g)["N_central"]
    assert math.isclose(N, 1.0 - g / math.sqrt(g * g + 1.0), rel_tol=1e-14)
    assert 0.0 < N < 1.0
    # GATE: B_centre = mu0 M (1 - N) must reproduce the INDEPENDENT cylinder
    # closed form (g = L/D, take a=1 -> D=2, L=2g)
    a = 1.0
    L = 2.0 * g * a            # D = 2a, so L/D = g
    centre = uniform_magnetized_cylinder_axial_field(M, a, L, 0.0)["centre"]
    assert math.isclose(centre, MU0 * M * (1.0 - N), rel_tol=1e-12)
    # long needle -> N -> 0
    assert cylinder_central_demag_factor(1e4)["N_central"] < 1e-6
    with pytest.raises(ValueError):
        cylinder_central_demag_factor(0.0)


def test_polarizability_sphere():
    mu_r = 5.0
    res = magnetic_polarizability_sphere(mu_r, A)
    beta_exact = 4.0 * math.pi * A ** 3 * (mu_r - 1.0) / (mu_r + 2.0)
    assert math.isclose(res["beta"], beta_exact, rel_tol=1e-14)
    assert math.isclose(res["M_avg_per_H0"],
                        3.0 * (mu_r - 1.0) / (mu_r + 2.0), rel_tol=1e-14)
    # GATE: mu_r = 1 (vacuum) -> beta = 0; mu_r -> inf -> M_avg/H0 -> 3
    assert magnetic_polarizability_sphere(1.0, A)["beta"] == 0.0
    huge = magnetic_polarizability_sphere(1e12, A)
    assert math.isclose(huge["M_avg_per_H0"], 3.0, rel_tol=1e-9)
    assert math.isclose(huge["beta"], 4.0 * math.pi * A ** 3, rel_tol=1e-9)
    with pytest.raises(ValueError):
        magnetic_polarizability_sphere(5.0, -A)


def test_coaxial_force():
    b_rem, L, gap = 1.2, 0.02, 0.01
    res = coaxial_magnets_axial_force(b_rem, A, L, gap)
    m = (b_rem / MU0) * math.pi * A * A * L
    s = gap + L
    assert math.isclose(res["moment"], m, rel_tol=1e-14)
    assert math.isclose(res["sep_s"], s, rel_tol=1e-14)
    assert math.isclose(res["force_N"], 3.0 * MU0 * m * m / (2.0 * math.pi * s ** 4),
                        rel_tol=1e-14)
    # GATE: dipole-dipole force scales as 1/s^4 -- doubling s drops F by 16x.
    # Compare two gaps giving s and 2s (independent of the prefactor).
    s1 = gap + L
    gap2 = 2.0 * s1 - L                     # so that gap2 + L = 2 s1
    res2 = coaxial_magnets_axial_force(b_rem, A, L, gap2)
    assert math.isclose(res["force_N"] / res2["force_N"], 16.0, rel_tol=1e-12)
    with pytest.raises(ValueError):
        coaxial_magnets_axial_force(b_rem, A, L, gap, mu_r=0.5)


def test_continuous_loop_slot62_coaxial_force_gap_sweep_contract():
    b_rem, radius, length = 1.2, 0.01, 0.02
    summary = coaxial_magnets_force_gap_sweep_summary(
        b_rem,
        radius,
        length,
        gaps=[0.01, 0.04, 0.10],
    )

    assert summary["status"] == "ok"
    assert summary["schema"] == "radia-ngsolve.coaxial-magnets-force-gap-sweep.v1"
    assert summary["checks"]["force_decreases_with_gap"] is True
    assert summary["checks"]["force_s4_invariant_ok"] is True
    assert summary["max_force_s4_invariant_rel_error"] < 1.0e-14
    assert [row["sep_s"] for row in summary["rows"]] == pytest.approx([0.03, 0.06, 0.12])
    assert summary["rows"][0]["force_N"] / summary["rows"][1]["force_N"] == pytest.approx(16.0)
    assert summary["rows"][1]["force_N"] / summary["rows"][2]["force_N"] == pytest.approx(16.0)
    assert summary["force_ratio_first_last"] == pytest.approx(256.0)


def test_halbach_segmentation():
    # exact special value: N=8 -> sin(pi/4)/(pi/4) = 0.900316...
    f8 = halbach_segmentation_factor(8)["seg_factor"]
    assert math.isclose(f8, math.sin(math.pi / 4.0) / (math.pi / 4.0), rel_tol=1e-14)
    assert math.isclose(f8, 0.9003163161571062, rel_tol=1e-12)
    # exact N=4 -> sin(pi/2)/(pi/2) = 2/pi
    assert math.isclose(halbach_segmentation_factor(4)["seg_factor"],
                        2.0 / math.pi, rel_tol=1e-14)
    # GATE: monotone approach to the ideal continuous Halbach (K_seg -> 1)
    seq = [halbach_segmentation_factor(n)["seg_factor"] for n in (4, 8, 16, 64)]
    assert all(b > a for a, b in zip(seq, seq[1:]))
    assert math.isclose(halbach_segmentation_factor(100000)["seg_factor"], 1.0,
                        rel_tol=1e-6)
    with pytest.raises(ValueError):
        halbach_segmentation_factor(2)


if __name__ == "__main__":
    test_sphere_field()
    test_cylinder_axial_field()
    test_annulus_axial_field()
    test_prism_demag_factors()
    test_cylinder_central_demag()
    test_polarizability_sphere()
    test_coaxial_force()
    test_halbach_segmentation()
    print("[OK] magnetized-body closed forms verified.")
