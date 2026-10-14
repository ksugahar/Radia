# -*- coding: utf-8 -*-
"""ngsolve.bem integral-equation (BEM) capability -- the OPEN integral solver, the radia-ngsolve
counterpart of a commercial high-frequency integral-equation (MoM) solver.

Verifies the Laplace SINGLE-LAYER BEM (the kappa=0 static limit of the high-frequency integral solver)
against the exact analytic capacitance of an isolated conducting sphere, C = 4 pi eps0 R (eps0=1 here).
"""
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

bem = pytest.importorskip("ngsolve.bem")
from radia_mcp.radia_ngsolve.bem_integral import (sphere_capacitance, single_layer_capacitance,
                                                  spheroid_capacitance, spheroid_capacitance_analytic,
                                                  helmholtz_exterior_point_source,
                                                  helmholtz_soundsoft_sphere_scattering,
                                                  mie_soundsoft_scattered,
                                                  helmholtz_cfie_soundsoft_sphere,
                                                  helmholtz_irregular_frequency_condition,
                                                  helmholtz_soundsoft_far_field,
                                                  mie_soundsoft_farfield, mie_soundsoft_cross_section,
                                                  helmholtz_soundhard_far_field,
                                                  mie_soundhard_farfield, mie_soundhard_cross_section,
                                                  helmholtz_impedance_far_field, mie_impedance_farfield,
                                                  helmholtz_sphere_radiation_impedance, sphere_radiation_ratio,
                                                  sphere_polarizability, conductor_polarizability,
                                                  single_layer_spectrum)


def test_spheroid_capacitance_prolate_oblate():
    # conducting SPHEROID capacitance via the Laplace single-layer BEM -- the non-spherical generalisation
    # of C = 4 pi eps0 R, exact-gated by the closed form (prolate C=4 pi cf/artanh(cf/a), oblate
    # C=4 pi cf/arcsin(cf/b), sphere 4 pi a; eps0=1). Strengthens the single-layer solver to a new geometry.
    import math
    # the analytic closed form: sphere limit, prolate, oblate disk limit
    assert abs(spheroid_capacitance_analytic(1.0, 1.0) - 4.0 * math.pi) < 1e-12        # sphere = 4 pi a
    assert abs(spheroid_capacitance_analytic(1.0, 0.6) - 9.150723) < 1e-4              # prolate
    assert abs(spheroid_capacitance_analytic(0.6, 1.0) - 10.841311) < 1e-4            # oblate
    assert abs(spheroid_capacitance_analytic(1e-4, 1.0) - 8.0) < 5e-3                  # oblate disk -> C=8 b
    with pytest.raises(ValueError):
        spheroid_capacitance_analytic(0.0, 1.0)
    # PROLATE 2:1.2 (a=1, b=0.6): the BEM single-layer reproduces the closed form
    pr = spheroid_capacitance(1.0, 0.6, maxh=0.18, order=3)
    assert pr["kind"] == "prolate"
    assert pr["rel_err"] < 5e-3, pr                                                    # curved order-3 BEM
    assert abs(pr["capacitance"] - spheroid_capacitance_analytic(1.0, 0.6)) < 5e-2
    # OBLATE (a=0.6, b=1.0): same solver, oblate closed form
    ob = spheroid_capacitance(0.6, 1.0, maxh=0.18, order=3)
    assert ob["kind"] == "oblate" and ob["rel_err"] < 5e-3, ob
    # the spheroid is MORE capacitive than the equal-volume... (qualitative: prolate C > inscribed sphere)
    assert pr["capacitance"] > 4.0 * math.pi * 0.6                                     # > sphere of radius b


def test_sphere_capacitance_matches_4pi():
    # unit conducting sphere -> C = 4 pi (eps0=1) via the Laplace single-layer indirect BEM
    res = sphere_capacitance(R=1.0, maxh=0.4, order=3)
    assert res["rel_err"] < 2e-3, res                # curved order-3 surface: ~machine precision
    assert abs(res["capacitance"] - res["analytic"]) < 1e-2


def test_capacitance_mesh_convergence():
    # the single-layer BEM capacitance converges to 4 pi as the curved surface is refined
    errs = [sphere_capacitance(R=1.0, maxh=h, order=3)["rel_err"] for h in (0.6, 0.4)]
    assert all(e < 5e-3 for e in errs)
    assert max(errs) < 1e-2


def test_capacitance_scales_with_radius():
    # C = 4 pi eps0 R is linear in R (geometric); a radius-2 sphere has twice the capacitance
    c1 = sphere_capacitance(R=1.0, maxh=0.4, order=3)["capacitance"]
    c2 = sphere_capacitance(R=2.0, maxh=0.8, order=3)["capacitance"]
    assert abs(c2 / c1 - 2.0) < 5e-3


def test_helmholtz_high_frequency_point_source():
    # HIGH-FREQUENCY Helmholtz single-layer BEM (the kappa>0 generalisation of the static capacitance):
    # an interior point source -> exterior field reproduced to ~machine precision (manufactured solution),
    # across ka = 1, 2, 4 (high-frequency).  The open integral solver vs the exact point-source field.
    for kappa in (1.0, 2.0, 4.0):
        res = helmholtz_exterior_point_source(kappa, a=1.0, maxh=0.5, order=3)
        assert res["rel_err"] < 1e-3, (kappa, res["rel_err"])    # measured ~1e-6
        assert abs(abs(res["u_exact"]) - 1.0 / (4 * 3.141592653589793 * 3.0)) < 1e-9  # |G| = 1/(4 pi r)


def test_helmholtz_soundsoft_mie_scattering():
    # HIGH-FREQUENCY scattering: an incident plane wave exp(i kappa z) on a SOUND-SOFT (Dirichlet) sphere,
    # solved by the ngsolve.bem single-layer BEM, vs the EXACT scalar Mie series -- a TRUE scattering
    # problem (real obstacle), not a manufactured solution.  ka < pi (single-layer irregular freq at pi).
    for kappa in (0.5, 1.0, 2.0):
        res = helmholtz_soundsoft_sphere_scattering(kappa, a=1.0, maxh=0.35, order=3, rext=2.0)
        assert abs(res["ka"] - kappa) < 1e-12
        for label in ("axial", "equator"):
            assert res[label]["rel_err"] < 1e-4, (kappa, label, res[label]["rel_err"])  # measured ~3e-7
            assert abs(res[label]["u_bem"] - res[label]["u_mie"]) < 1e-4


def test_mie_soundsoft_low_ka_limit():
    # Rayleigh (ka->0) limit: the n=0 monopole dominates, scattered field ~ -(j_0(ka)/h_0(ka)) h_0(k r);
    # the Mie helper is finite and the BEM reproduces it (small-ka monopole scattering).
    u = mie_soundsoft_scattered(0.3, 1.0, 0.0, 0.0, 3.0)
    assert abs(u) > 0.0 and abs(u.imag) > 0.0           # complex, non-trivial scattered field
    res = helmholtz_soundsoft_sphere_scattering(0.3, a=1.0, maxh=0.4, order=3, rext=3.0)
    assert res["axial"]["rel_err"] < 1e-3


def test_cfie_scattering_across_irregular_frequency():
    # COMBINED-FIELD (Brakhage-Werner) sound-soft scattering is robust ACROSS the single-layer's irregular
    # frequency ka=pi (the helmholtz_soundsoft_sphere_scattering single-layer is limited to ka<pi).  The
    # CFIE reproduces the scalar Mie series to ~machine precision for all ka including ka=pi.
    import math
    for kappa in (2.0, math.pi, 4.0):                   # 4.0 and pi are >= the single-layer irregular freq
        res = helmholtz_cfie_soundsoft_sphere(kappa, a=1.0, maxh=0.35, order=3, rext=2.0)
        assert abs(res["ka"] - kappa) < 1e-12
        for label in ("axial", "equator"):
            assert res[label]["rel_err"] < 1e-4, (kappa, label, res[label]["rel_err"])  # measured ~1e-6


def test_combined_field_well_conditioned_at_irregular_frequency():
    # The single-layer V is RANK-DEFICIENT at the interior Dirichlet eigenvalue ka=pi (j_0(pi)=0) -> cond(V)
    # spikes; the combined-field operator (1/2 M + CFO) is uniquely solvable for all kappa -> bounded cond.
    import math
    off = helmholtz_irregular_frequency_condition(2.5, a=1.0, maxh=0.9, order=0)
    res = helmholtz_irregular_frequency_condition(math.pi, a=1.0, maxh=0.9, order=0)
    # single-layer condition number SPIKES at ka=pi (measured ~4.6e3 vs ~15 off-resonance)
    assert res["cond_single_layer"] > 20.0 * off["cond_single_layer"], \
        (off["cond_single_layer"], res["cond_single_layer"])
    # combined-field stays well-conditioned and FLAT across the irregular frequency (measured ~3.35)
    assert res["cond_combined_field"] < 10.0
    assert abs(res["cond_combined_field"] - off["cond_combined_field"]) < 1.0
    assert res["cond_combined_field"] < res["cond_single_layer"]


def test_soundsoft_far_field_and_optical_theorem():
    # The RCS apparatus: the FAR-FIELD scattering amplitude f(theta) from the combined-field density matches
    # the scalar Mie far-field to ~machine precision, AND the OPTICAL THEOREM sigma_t=(4pi/k)Im f(0)==sigma_sca.
    res = helmholtz_soundsoft_far_field(2.0, a=1.0, maxh=0.3, order=3)
    for fb, fm, e in zip(res["f_bem"], res["f_mie"], res["rel_err"]):
        assert e < 1e-4, (fb, fm, e)                       # measured ~1e-7 across theta=0..180
    # optical theorem: forward Im f(0) gives the total cross-section == Mie sigma_sca
    assert res["optical_rel_err"] < 1e-4, res              # measured ~2.5e-7
    assert res["sigma_optical"] > 0.0 and abs(res["sigma_optical"] - res["sigma_mie"]) < 1e-3


def test_soundhard_far_field_and_optical_theorem():
    # #225: the RIGID-body RCS apparatus -- the SOUND-HARD (Neumann) far-field scattering amplitude f(theta)
    # from the direct (K-1/2 M)u=V sigma solve matches the rigid Mie far-field, AND the optical theorem
    # sigma_t=(4pi/k)Im f(0)==sigma_sca. The Neumann DUAL of test_soundsoft_far_field_and_optical_theorem;
    # the open ngsolve.bem RCS of a rigid sphere (the open counterpart of a commercial high-frequency MoM RCS).
    res = helmholtz_soundhard_far_field(2.0, a=1.0, maxh=0.3, order=3)
    for fb, fm, e in zip(res["f_bem"], res["f_mie"], res["rel_err"]):
        assert e < 1e-4, (fb, fm, e)                       # measured ~1-7e-7 across theta=0..180
    assert res["optical_rel_err"] < 1e-4, res              # measured ~3e-7
    assert res["sigma_optical"] > 0.0 and abs(res["sigma_optical"] - res["sigma_mie"]) < 1e-3
    # rigid Mie analytic gates: the forward far-field + cross-section satisfy the optical theorem exactly
    import math as _m
    for kappa in (1.0, 2.0, 4.0):
        f0 = mie_soundhard_farfield(kappa, 1.0, 0.0)
        sig = mie_soundhard_cross_section(kappa, 1.0)
        assert abs(4 * _m.pi / kappa * f0.imag - sig) < 1e-9 * sig    # optical theorem (analytic identity)
        assert sig > 0.0
        assert abs(mie_soundhard_farfield(kappa, 1.0, _m.pi)) > 0.0   # backscatter finite


def test_impedance_far_field():
    # #230: the IMPEDANCE (Robin) sphere scattering du/dn + i k beta u = 0 -- the INTERPOLATION between the
    # sound-soft (beta->inf) and sound-hard (beta=0) limits, via the SAME direct (K-1/2 M + i k beta V)u=V g
    # solve reusing the sound-hard V/K/M operators. The impedance / coated-obstacle member of the bem family.
    import math as _m
    # ANALYTIC limits: the impedance Mie ratio s_n -> j_n'/h_n' (sound-hard) at beta=0, -> j_n/h_n (soft) at beta->inf
    for ka in (1.0, 2.0):
        assert abs(mie_impedance_farfield(ka, 1.0, 0.0, 0.7) - mie_soundhard_farfield(ka, 1.0, 0.7)) < 1e-12
        assert abs(mie_impedance_farfield(ka, 1.0, 1e7, 0.7) - mie_soundsoft_farfield(ka, 1.0, 0.7)) < 1e-4
    # the BEM Robin solve reproduces the impedance Mie far-field
    res = helmholtz_impedance_far_field(2.0, a=1.0, beta=0.5, maxh=0.3, order=3)
    assert abs(res["beta"] - 0.5) < 1e-12
    for fb, fm, e in zip(res["f_bem"], res["f_mie"], res["rel_err"]):
        assert e < 1e-4, (fb, fm, e)                       # measured ~1-8e-7 across theta=0..180


def test_helmholtz_sphere_radiation_impedance():
    # #220: the sphere RADIATION impedance (velocity-driven Neumann) -- the RADIATION complement of the
    # scattering family, reusing the sound-hard direct Calderon (K - 1/2 M) u = V sigma. The surface ratio
    # u/sigma = h_n(ka)/(k h_n'(ka)) for the monopole (n=0) and dipole (n=1), reproduced to ~1e-7.
    import cmath
    # the analytic monopole ratio = a/(ika-1); the radiation resistance R/(rho c)=(ka)^2/(1+(ka)^2)
    for ka in (0.5, 1.0, 2.0):
        an = sphere_radiation_ratio(ka, 1.0, 0)
        assert abs(an - 1.0 / (1j * ka - 1.0)) < 1e-12
        # z = j omega rho (u/sigma) -> for the unit-rho-c normalisation the monopole resistance is
        # (ka)^2/(1+(ka)^2): Re of (j ka)*(a/(ika-1))/a ... checked via the known closed form
        znorm = 1j * ka * an / 1.0                                   # j ka * (u/sigma), a=1
        assert abs(znorm.real - ka ** 2 / (1 + ka ** 2)) < 1e-9      # radiation resistance R/(rho c)
    # the BEM velocity-driven solve reproduces the radiation ratio for both modes
    for kappa in (0.5, 1.0, 2.0):
        r = helmholtz_sphere_radiation_impedance(kappa, a=1.0, maxh=0.3, order=3)
        assert abs(r["ka"] - kappa) < 1e-12
        for mode in ("monopole", "dipole"):
            assert r[mode]["rel_err"] < 1e-4, (kappa, mode, r[mode]["rel_err"])   # measured ~1e-7
            assert abs(r[mode]["ratio_bem"] - r[mode]["ratio_analytic"]) < 1e-4


def test_mie_farfield_and_cross_section_consistency():
    # standalone analytic gates: the forward far-field and the cross-section satisfy the optical theorem
    # exactly (sigma_sca = (4pi/k) Im f(0) for the lossless sound-soft sphere) across ka.
    import math as _m
    for kappa in (1.0, 2.0, 4.0):
        f0 = mie_soundsoft_farfield(kappa, 1.0, 0.0)
        sig = mie_soundsoft_cross_section(kappa, 1.0)
        assert abs(4 * _m.pi / kappa * f0.imag - sig) < 1e-9 * sig    # optical theorem (analytic identity)
        assert sig > 0.0
        # back-scatter (theta=pi) is finite and nonzero
        assert abs(mie_soundsoft_farfield(kappa, 1.0, _m.pi)) > 0.0


def test_single_layer_operator_spectrum():
    # eigenvalues of the Laplace single-layer operator on the unit sphere: exact V Y_l = Y_l/(2l+1),
    # so the spectrum is 1, 1/3 (x3), 1/5 (x5) descending -- the operator companion of the capacitance solve.
    res = single_layer_spectrum(n_eig=9, maxh=0.6, order=1)
    ev = res["eigenvalues"]
    assert ev[0] == pytest.approx(1.0, rel=5e-3)               # l=0 monopole (x1)
    for v in ev[1:4]:
        assert v == pytest.approx(1.0 / 3.0, rel=5e-3)         # l=1 (x3)
    for v in ev[4:9]:
        assert v == pytest.approx(1.0 / 5.0, rel=5e-3)         # l=2 (x5)
    assert all(a >= b - 1e-9 for a, b in zip(ev, ev[1:]))      # descending


def test_sphere_polarizability_matches_4piR3():
    # conducting sphere in a UNIFORM applied field -> induced-dipole polarizability alpha = 4 pi eps0 R^3
    # (eps0=1), via the single-layer BEM with a COORDINATE datum (V sigma = z) instead of the constant-1
    # capacitance datum. The applied-field RESPONSE complement of the fixed-potential self-capacitance.
    import math
    res = sphere_polarizability(R=1.0, maxh=0.3, order=3)
    assert abs(res["analytic"] - 4.0 * math.pi) < 1e-12          # 4 pi R^3, R=1
    assert res["rel_err"] < 2e-3, res                            # curved order-3 BEM, measured ~machine
    assert abs(res["qnet"]) < 1e-5                               # neutral conductor (zero net charge)


def test_polarizability_scales_with_volume():
    # alpha = 4 pi eps0 R^3 is CUBIC in R (a volume-like dipole response), unlike capacitance C ~ R.
    a1 = sphere_polarizability(R=1.0, maxh=0.3, order=3)["alpha"]
    a2 = sphere_polarizability(R=2.0, maxh=0.6, order=3)["alpha"]
    assert abs(a2 / a1 - 8.0) < 2e-2                             # (2)^3 = 8


def test_polarizability_isotropic_x_vs_z():
    # the sphere is isotropic: the x-axis and z-axis applied-field responses are equal.
    import math
    mesh_kw = dict(maxh=0.3, order=3)
    az = sphere_polarizability(R=1.0, **mesh_kw)["alpha"]
    # rebuild the same mesh and drive along x -> same polarizability
    import ngsolve as ng
    from netgen.occ import Sphere, Pnt, OCCGeometry
    mesh = ng.Mesh(OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0)).GenerateMesh(maxh=0.3)).Curve(4)
    ax = conductor_polarizability(mesh, axis="x", order=3)["alpha"]
    assert abs(ax - az) / az < 5e-3
    assert abs(az - 4.0 * math.pi) / (4.0 * math.pi) < 2e-3
