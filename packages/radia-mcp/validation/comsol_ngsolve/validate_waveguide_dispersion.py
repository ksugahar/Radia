r"""Waveguide dispersion -- propagation above cutoff + evanescence below (#65) -- test.

beta=(2pi/c)sqrt(f^2-fc^2); lambda_g=lambda0/sqrt(1-(fc/f)^2); v_p=c/sqrt(1-(fc/f)^2);
v_g=c sqrt(1-(fc/f)^2); v_p*v_g=c^2; below cutoff alpha=(2pi/c)sqrt(fc^2-f^2). Closed-form
identities + a numerical d omega/d beta check + an FE-anchored cutoff (the propagation sequel to
the #53/#62 cutoff eigenvalues)."""
import cmath
import math
import os
import sys

import pytest

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

from radia_mcp.radia_ngsolve.waveguide import (rectangular_waveguide_cutoff, cutoff_frequency,
                                               waveguide_dispersion, guide_wavelength,
                                               waveguide_evanescent_attenuation,
                                               waveguide_dielectric_slab_sparams,
                                               waveguide_cascade_sparams,
                                               waveguide_offset_short_s11, C0)


def validate_dielectric_slab_sparams():
    # #216: the dielectric-slab-loaded TE10 waveguide 2-port S-matrix (transmission-line) -- the MISMATCHED
    # sequel to the matched line (waveguide_dispersion / the matched-port #131). Exact gates: lossless
    # unitarity, slab->0 / eps_r=1 matched limits, quarter-wave max reflection, the air-region beta.
    a = 0.02286                                          # WR-90 width; TE10 fc = c/2a = 6.557 GHz
    s = waveguide_dielectric_slab_sparams(10e9, a, 2.2, 0.010)
    # LOSSLESS UNITARITY |S11|^2 + |S21|^2 = 1 (exact)
    assert math.isclose(s["unitarity"], 1.0, rel_tol=1e-12)
    # fc + propagation constants: beta_slab > beta_air (denser medium), Gamma < 0 for eps_r > 1
    assert math.isclose(s["fc"], 0.5 * C0 / a, rel_tol=1e-12)
    assert s["beta_slab"] > s["beta_air"] > 0.0 and s["gamma"] < 0.0
    # MATCHED limits: a zero-length slab and eps_r=1 both give S11=0, |S21|=1
    for sm in (waveguide_dielectric_slab_sparams(10e9, a, 2.2, 0.0),
               waveguide_dielectric_slab_sparams(10e9, a, 1.0, 0.010)):
        assert abs(sm["S11"]) < 1e-12 and math.isclose(sm["S21_mag"], 1.0, rel_tol=1e-12)
    # QUARTER-WAVE slab (theta=pi/2) maximises reflection at |S11| = 2|Gamma|/(1+Gamma^2)
    g = s["gamma"]
    dq = (math.pi / 2.0) / s["beta_slab"]
    sq = waveguide_dielectric_slab_sparams(10e9, a, 2.2, dq)
    assert math.isclose(sq["theta"], math.pi / 2.0, rel_tol=1e-9)
    assert math.isclose(sq["S11_mag"], 2.0 * abs(g) / (1.0 + g * g), rel_tol=1e-9)
    # the air-region beta equals the matched-line waveguide_dispersion beta (consistency with #65/#131)
    assert math.isclose(s["beta_air"], waveguide_dispersion(10e9, s["fc"])["beta"], rel_tol=1e-12)
    # validation
    with pytest.raises(ValueError):
        waveguide_dielectric_slab_sparams(5e9, a, 2.2, 0.010)        # below TE10 cutoff
    with pytest.raises(ValueError):
        waveguide_dielectric_slab_sparams(10e9, a, 0.0, 0.010)       # eps_r must be > 0


def validate_cascade_and_offset_short():
    # #221: the N-section waveguide ABCD CASCADE (generalises the single slab #216) + the OFFSET-SHORT
    # 1-port. Exact gates: one-section cascade == the slab to machine precision; det(ABCD)=1
    # (reciprocity); lossless unitarity; matched air line; transparent half-wave slab; the offset
    # short == a short-terminated air line == -exp(-2j beta_air d). (LIVE-anchored: the offset-short
    # S11(f) was reproduced by a full-wave 1-port port+PEC solve to ~3e-5, |S11|-1 ~ 5e-15.)
    a = 0.02286                                              # WR-90; TE10 fc = c/2a = 6.557 GHz
    f, er, d = 10e9, 2.2, 0.010
    # (1) ONE section reproduces waveguide_dielectric_slab_sparams to machine precision
    slab = waveguide_dielectric_slab_sparams(f, a, er, d)
    cas1 = waveguide_cascade_sparams(f, a, [(d, er)])
    assert abs(cas1["S11"] - slab["S11"]) < 1e-12 and abs(cas1["S21"] - slab["S21"]) < 1e-12
    assert cas1["n_sections"] == 1
    # reciprocity det(ABCD)=1 and lossless unitarity
    assert math.isclose(cas1["abcd_det"].real, 1.0, rel_tol=1e-12) and abs(cas1["abcd_det"].imag) < 1e-12
    assert math.isclose(cas1["unitarity"], 1.0, rel_tol=1e-12)
    # (2) a MATCHED air section (eps_r=1) -> S11=0, |S21|=1
    air = waveguide_cascade_sparams(f, a, [(d, 1.0)])
    assert abs(air["S11"]) < 1e-12 and math.isclose(air["S21_mag"], 1.0, rel_tol=1e-12)
    # (3) a HALF-WAVE slab (theta=pi) is transparent: S11=0, |S21|=1
    dhalf = math.pi / slab["beta_slab"]
    half = waveguide_cascade_sparams(f, a, [(dhalf, er)])
    assert abs(half["S11"]) < 1e-9 and math.isclose(half["S21_mag"], 1.0, rel_tol=1e-9)
    # (4) MULTI-section air|slab|air: lossless-unitary, reciprocal, and (lossless) |S11|,|S21| match
    #     the bare slab -- the air sections are reference-plane only
    multi = waveguide_cascade_sparams(f, a, [(0.005, 1.0), (d, er), (0.005, 1.0)])
    assert multi["n_sections"] == 3 and math.isclose(multi["unitarity"], 1.0, rel_tol=1e-12)
    assert math.isclose(multi["abcd_det"].real, 1.0, rel_tol=1e-9) and abs(multi["abcd_det"].imag) < 1e-9
    assert math.isclose(multi["S11_mag"], slab["S11_mag"], rel_tol=1e-9)
    assert math.isclose(multi["S21_mag"], slab["S21_mag"], rel_tol=1e-9)
    # (5) OFFSET SHORT: |S11|=1, S11 = -exp(-2j beta_air d)
    os_ = waveguide_offset_short_s11(f, a, d)
    assert math.isclose(os_["S11_mag"], 1.0, rel_tol=1e-12)
    assert abs(os_["S11"] - (-cmath.exp(-2j * os_["beta_air"] * d))) < 1e-12
    # a short-terminated air line (Gamma_L=-1) reproduces it: Gamma_in = -S21^2 (S11=S22=0)
    gin = air["S11"] + air["S21"] * air["S21"] * (-1.0) / (1.0 - air["S11"] * (-1.0))
    assert abs(gin - os_["S11"]) < 1e-9
    # the offset-short beta_air equals the matched-line dispersion beta (consistency with #65/#216)
    assert math.isclose(os_["beta_air"], waveguide_dispersion(f, os_["fc"])["beta"], rel_tol=1e-12)
    # validation
    with pytest.raises(ValueError):
        waveguide_cascade_sparams(5e9, a, [(d, er)])            # below TE10 cutoff
    with pytest.raises(ValueError):
        waveguide_cascade_sparams(f, a, [])                     # empty cascade
    with pytest.raises(ValueError):
        waveguide_offset_short_s11(5e9, a, d)                   # below cutoff


def validate_dispersion_identities():
    fc = rectangular_waveguide_cutoff(0.02286, 0.01016, 1, 0)     # WR-90 TE10
    for f in (8e9, 10e9, 12e9, 15e9):
        d = waveguide_dispersion(f, fc)
        # the reciprocal velocity identity v_p * v_g = c^2 (exact)
        assert math.isclose(d["v_phase"] * d["v_group"], C0 ** 2, rel_tol=1e-12)
        # phase super-luminal, group sub-luminal
        assert d["v_phase"] > C0 > d["v_group"] > 0.0
        # guide-wavelength relation 1/lg^2 = 1/l0^2 - 1/lc^2 (exact) and lg > l0
        lam0, lam_c = C0 / f, C0 / fc
        assert math.isclose(1.0 / d["lambda_g"] ** 2, 1.0 / lam0 ** 2 - 1.0 / lam_c ** 2, rel_tol=1e-12)
        assert d["lambda_g"] > lam0
        assert math.isclose(guide_wavelength(f, fc), d["lambda_g"], rel_tol=1e-12)
        # group velocity == d omega/d beta from omega(beta) = c sqrt(beta^2 + kc^2) (central diff)
        kc = 2.0 * math.pi * fc / C0
        db = d["beta"] * 1e-6
        w = lambda bb: C0 * math.sqrt(bb * bb + kc * kc)
        vg_num = (w(d["beta"] + db) - w(d["beta"] - db)) / (2.0 * db)
        assert math.isclose(vg_num, d["v_group"], rel_tol=1e-6)


def validate_cutoff_and_high_frequency_limits():
    fc = 6.5571e9
    # propagation branch raises at/below cutoff; evanescent branch raises at/above it
    with pytest.raises(ValueError):
        waveguide_dispersion(fc, fc)
    with pytest.raises(ValueError):
        waveguide_dispersion(0.9 * fc, fc)
    with pytest.raises(ValueError):
        waveguide_evanescent_attenuation(fc, fc)
    # just above cutoff: v_g -> 0, v_p -> inf, beta -> 0
    d = waveguide_dispersion(1.0001 * fc, fc)
    assert d["v_group"] < 0.02 * C0 and d["v_phase"] > 50.0 * C0
    assert waveguide_dispersion(1.000001 * fc, fc)["beta"] < 1.0
    # high frequency: both velocities -> c
    d = waveguide_dispersion(100.0 * fc, fc)
    assert math.isclose(d["v_phase"], C0, rel_tol=1e-3)
    assert math.isclose(d["v_group"], C0, rel_tol=1e-3)


def validate_evanescent_limits():
    fc = 6.5571e9
    kc = 2.0 * math.pi * fc / C0
    # alpha -> k_c as f -> 0, and alpha -> 0 as f -> fc
    assert math.isclose(waveguide_evanescent_attenuation(1e-3 * fc, fc), kc, rel_tol=1e-3)
    assert waveguide_evanescent_attenuation(0.9999 * fc, fc) < 0.02 * kc
    assert waveguide_evanescent_attenuation(0.999999 * fc, fc) < 1.0     # continuity at cutoff


def validate_dispersion_fe_anchor():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen")
    from ngsolve import Mesh, TaskManager
    from netgen.occ import OCCGeometry, WorkPlane
    from radia_mcp.radia_ngsolve.waveguide import helmholtz_cutoff_wavenumbers_2d

    a, b = 0.02286, 0.01016
    rect = WorkPlane().Rectangle(a, b).Face(); rect.edges.name = "wall"
    mesh = Mesh(OCCGeometry(rect, dim=2).GenerateMesh(maxh=min(a, b) / 12))
    with TaskManager():
        kc_fe = helmholtz_cutoff_wavenumbers_2d(mesh, 1, bc="neumann")[0]

    fc_fe = cutoff_frequency(kc_fe)
    fc_an = rectangular_waveguide_cutoff(a, b, 1, 0)
    assert math.isclose(fc_fe, fc_an, rel_tol=5e-3)                 # FE cutoff matches analytic
    # the dispersion built on the FE cutoff agrees with the analytic-cutoff dispersion
    d_fe = waveguide_dispersion(10e9, fc_fe)
    d_an = waveguide_dispersion(10e9, fc_an)
    assert math.isclose(d_fe["v_group"], d_an["v_group"], rel_tol=1e-2)
    assert math.isclose(d_fe["lambda_g"], d_an["lambda_g"], rel_tol=1e-2)


if __name__ == "__main__":
    validate_dielectric_slab_sparams()
    validate_cascade_and_offset_short()
    validate_dispersion_identities()
    validate_cutoff_and_high_frequency_limits()
    validate_evanescent_limits()
    validate_dispersion_fe_anchor()
    print("[OK] waveguide dispersion: slab/cascade/offset-short S-params, v_p*v_g=c^2, FE-anchored fc.")
